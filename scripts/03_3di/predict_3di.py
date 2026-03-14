"""
predict_3di.py - Predict 3Di sequences from amino acid FASTA using ProstT5.

Uses AutoModelForSeq2SeqLM with model.generate() for AA -> 3Di translation.
This is the officially supported full-model inference mode from the ProstT5
HuggingFace page (Rostlab/ProstT5), requiring no external CNN weights.

ProstT5 is a T5 encoder-decoder model fine-tuned on 17M AlphaFoldDB structures.
Direction of translation is controlled by a prefix token:
    "<AA2fold>"  : amino acid sequence -> 3Di tokens  (used here)
    "<fold2AA>"  : 3Di tokens -> amino acid sequence  (inverse folding)

AA sequences must be UPPER-CASE. 3Di output is LOWER-CASE (Foldseek convention).
Rare/ambiguous amino acids (U, Z, O, B) are replaced with X before tokenisation.

Dependencies:
    pip install "torch>=2.6" --index-url https://download.pytorch.org/whl/cu124
    pip install transformers sentencepiece

    Requires CUDA driver >= 12.4 (Phoenix: --gres=gg:g4:1).
    See scripts/03_3di/README for full installation notes.

Environment:
    export HF_HOME="$PROJECT_ROOT/.hf_cache"   # cache model outside conda dir
    export HF_HUB_DISABLE_XET=1                # disable unstable xet protocol

Usage:
    # Both classes (default paths from config):
    python scripts/03_3di/predict_3di.py

    # Explicit paths:
    python scripts/03_3di/predict_3di.py \\
        --class1 data/processed/class1_sequences.fasta \\
        --class2 data/processed/class2_sequences.fasta

    # Skip one class:
    python scripts/03_3di/predict_3di.py --no-class2

    # Reduce batch size if OOM:
    python scripts/03_3di/predict_3di.py --batch-size 8

Output:
    data/processed/class1_3di.fasta
    data/processed/class2_3di.fasta

    Headers are preserved from input FASTA.
    Sequences are lowercase 3Di tokens (a-y alphabet, 20 tokens).
    Sequences exceeding 1000 aa are skipped (ProstT5 hard limit).

Notes:
    - Model (~11GB) is downloaded on first run and cached at $HF_HOME.
      Subsequent runs load from cache — no re-download needed.
    - Checkpoint written after each batch — safe to kill and resume.
      Checkpoint is removed automatically on successful completion.
    - Batching respects both --batch-size and max_residues=4000 (sum of
      sequence lengths per batch) to avoid GPU OOM. Sequences are sorted
      longest-first to minimise padding waste.
"""

import sys
import re
import json
import argparse
import logging
from pathlib import Path

import torch

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

MODEL_ID     = "Rostlab/ProstT5"
MAX_SEQ_LEN  = 1000   # hard limit from ProstT5 paper
MAX_RESIDUES = 4000   # max total residues per batch to avoid OOM


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Return list of (header, sequence) from a FASTA file."""
    records = []
    header, seq_parts = None, []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())

    if header is not None:
        records.append((header, "".join(seq_parts)))

    log.info(f"  Parsed {len(records):,} sequences from {path.name}")
    return records


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    """Write (header, sequence) pairs to a FASTA file."""
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
    log.info(f"  Written {len(records):,} sequences -> {path}")


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------

def load_model(device: torch.device):
    """
    Load ProstT5 tokenizer and full encoder-decoder model.
    Uses fp16 on GPU, fp32 on CPU.
    """
    from transformers import T5Tokenizer, AutoModelForSeq2SeqLM

    log.info(f"Loading ProstT5 tokenizer ({MODEL_ID})...")
    tokenizer = T5Tokenizer.from_pretrained(MODEL_ID, do_lower_case=False)

    log.info(f"Loading ProstT5 model ({MODEL_ID})...")
    model = AutoModelForSeq2SeqLM.from_pretrained(MODEL_ID).to(device)

    if device.type == "cuda":
        model.half()   # fp16 on GPU as recommended by ProstT5 docs
        log.info("  Model loaded in fp16 on GPU")
    else:
        model.full()
        log.warning("  Model loaded in fp32 on CPU — this will be very slow")

    model.eval()
    return tokenizer, model


# ---------------------------------------------------------------------------
# 3Di prediction
# ---------------------------------------------------------------------------

def make_batches(records: list[tuple[str, str]],
                 batch_size: int) -> list[list[tuple[str, str]]]:
    """
    Split records into batches respecting both batch_size and MAX_RESIDUES.
    Sorts by length (longest first) to minimise padding waste.
    """
    # Sort by descending length for efficient batching
    sorted_records = sorted(records, key=lambda x: len(x[1]), reverse=True)

    batches = []
    current_batch = []
    current_residues = 0

    for header, seq in sorted_records:
        seq_len = len(seq)
        if (len(current_batch) >= batch_size or
                current_residues + seq_len > MAX_RESIDUES) and current_batch:
            batches.append(current_batch)
            current_batch = []
            current_residues = 0
        current_batch.append((header, seq))
        current_residues += seq_len

    if current_batch:
        batches.append(current_batch)

    return batches


def predict_3di_batch(sequences: list[str],
                      tokenizer,
                      model,
                      device: torch.device) -> list[str]:
    """
    Predict 3Di tokens for a batch of AA sequences using generate().

    Prepends "<AA2fold>" prefix, replaces rare AAs with X, adds spaces
    between residues (required by ProstT5 tokenizer).

    Returns list of lowercase 3Di strings, one per input sequence.
    """
    # Replace rare/ambiguous AAs with X, add spaces, prepend direction token
    cleaned  = [re.sub(r"[UZOB]", "X", seq) for seq in sequences]
    prefixed = ["<AA2fold> " + " ".join(list(s)) for s in cleaned]

    # Tokenise
    ids = tokenizer(
        prefixed,
        add_special_tokens=True,
        padding="longest",
        return_tensors="pt",
    ).to(device)

    # Generate 3Di tokens (greedy decoding)
    with torch.no_grad():
        outputs = model.generate(
            input_ids=ids.input_ids,
            attention_mask=ids.attention_mask,
            max_new_tokens=max(len(s) for s in sequences) + 2,
            do_sample=False,        # greedy
        )

    results = []
    for i, seq in enumerate(sequences):
        # Decode output tokens, skip special tokens, remove spaces, lowercase
        decoded = tokenizer.decode(outputs[i], skip_special_tokens=True)
        tdi_seq = decoded.replace(" ", "").lower()

        results.append(tdi_seq)

    return results


def predict_3di(records: list[tuple[str, str]],
                tokenizer,
                model,
                device: torch.device,
                batch_size: int,
                checkpoint_path: Path) -> list[tuple[str, str]]:
    """
    Predict 3Di for all records with checkpointing for resumability.
    Skips sequences already in checkpoint.
    Writes checkpoint after each batch.
    Returns list of (header, 3di_sequence) in original order.
    """
    # Load existing checkpoint
    done: dict[str, str] = {}
    if checkpoint_path.exists():
        done = json.loads(checkpoint_path.read_text())
        log.info(f"  Resuming: {len(done):,} sequences already predicted")

    # Filter out already done and too-long sequences
    skipped = [h for h, s in records if len(s) > MAX_SEQ_LEN]
    if skipped:
        log.warning(
            f"  Skipping {len(skipped)} sequences exceeding {MAX_SEQ_LEN} aa: "
            + ", ".join(h[:50] for h in skipped[:3])
            + ("..." if len(skipped) > 3 else "")
        )

    todo = [
        (h, s) for h, s in records
        if h not in done and len(s) <= MAX_SEQ_LEN
    ]

    n_total = len(records) - len(skipped)
    log.info(f"  To predict: {len(todo):,} sequences (batch_size={batch_size}, "
             f"max_residues={MAX_RESIDUES})")

    batches = make_batches(todo, batch_size)
    log.info(f"  Split into {len(batches)} batches")

    for batch_idx, batch in enumerate(batches):
        headers = [h for h, _ in batch]
        seqs    = [s for _, s in batch]

        predictions = predict_3di_batch(seqs, tokenizer, model, device)

        for header, tdi in zip(headers, predictions):
            done[header] = tdi

        # Save checkpoint after every batch
        checkpoint_path.write_text(json.dumps(done))

        log.info(
            f"  Batch {batch_idx+1}/{len(batches)} done | "
            f"total: {len(done):,}/{n_total:,} sequences"
        )

    # Reconstruct in original input order
    results = []
    for header, seq in records:
        if header in done:
            results.append((header, done[header]))
        # Sequences that were skipped (too long) are omitted

    return results


# ---------------------------------------------------------------------------
# Per-file processing
# ---------------------------------------------------------------------------

def process_file(input_path: Path,
                 output_path: Path,
                 tokenizer,
                 model,
                 device: torch.device,
                 batch_size: int) -> None:
    """Predict 3Di for all sequences in input_path, write to output_path."""

    if output_path.exists():
        log.info(f"[SKIP] Output already exists: {output_path}")
        return

    checkpoint_path = output_path.with_suffix(".checkpoint.json")
    records = parse_fasta(input_path)

    results = predict_3di(
        records, tokenizer, model, device, batch_size, checkpoint_path
    )

    write_fasta(results, output_path)

    # Remove checkpoint on successful completion
    if checkpoint_path.exists():
        checkpoint_path.unlink()
        log.info("  Checkpoint removed")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Predict 3Di sequences from AA FASTA using ProstT5."
    )
    parser.add_argument(
        "--class1", type=Path,
        default=Path(PATHS["data_processed"]) / "class1_sequences.fasta",
    )
    parser.add_argument(
        "--class2", type=Path,
        default=Path(PATHS["data_processed"]) / "class2_sequences.fasta",
    )
    parser.add_argument("--no-class1", action="store_true", help="Skip Class I")
    parser.add_argument("--no-class2", action="store_true", help="Skip Class II")
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(PATHS["data_processed"]),
        help="Output directory (default: data/processed/)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=16,
        help="Max sequences per batch (also capped by max_residues=4000, default: 16)"
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Device: {device}")
    if device.type == "cpu":
        log.warning(
            "No GPU detected — prediction will be extremely slow. "
            "Request a GPU node with: #SBATCH --gres=gpu:1"
        )

    # Verify PyTorch sees the GPU
    if device.type == "cuda":
        log.info(f"  GPU: {torch.cuda.get_device_name(0)}")
        log.info(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    # Load model once, reuse for both classes
    tokenizer, model = load_model(device)

    to_run = []
    if not args.no_class1:
        to_run.append((args.class1, args.outdir / "class1_3di.fasta", "Class I"))
    if not args.no_class2:
        to_run.append((args.class2, args.outdir / "class2_3di.fasta", "Class II"))

    if not to_run:
        log.error("Nothing to run — both --no-class1 and --no-class2 specified.")
        sys.exit(1)

    for input_path, output_path, label in to_run:
        if not input_path.exists():
            log.error(f"{label}: input not found: {input_path}")
            log.error("Run scripts/01_search/extract_sequences.py first.")
            continue

        log.info(f"\n{'='*50}")
        log.info(f" {label}: {input_path.name} -> {output_path.name}")
        log.info(f"{'='*50}")
        process_file(input_path, output_path, tokenizer, model, device, args.batch_size)

    log.info("\nDone.")
    log.info("Next step — align 3Di FASTAs:")
    log.info("  bash scripts/02_align/run_mafft.sh \\")
    log.info("      data/processed/class1_3di.fasta \\")
    log.info("      results/alignments/class1_3di_aligned.fasta --threads 8")


if __name__ == "__main__":
    main()