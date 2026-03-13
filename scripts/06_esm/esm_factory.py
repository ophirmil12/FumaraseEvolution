"""
esm_factory.py - ESM-2 protein sequence embeddings.

Wraps ESMFactory with a CLI for batch embedding of FASTA files,
saving results as torch .pt dicts for downstream use.

Usage:
    python scripts/06_esm/esm_factory.py \\
        --fasta  data/processed/class1_sequences.fasta \\
        --output data/processed/embeddings/class1_embeddings.pt \\
        --embedding-size 1280 \\
        --batch-size 32

Output:
    .pt file containing {seq_id: mean_embedding_tensor} (Shape: D)
    Saved with torch.save(), loaded with torch.load(map_location="cpu")

Notes:
    - Sequences longer than 1022 aa are truncated (ESM-2 hard limit)
    - Checkpoint written after each batch — safe to kill and resume
    - Checkpoint removed automatically on successful completion
"""

import sys
import argparse
import logging
from pathlib import Path
import esm
from torch.utils.data import DataLoader, Dataset
from typing import List, Tuple, Dict, Iterator
from tqdm import tqdm

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

class ProteinDataset(Dataset):
    """Simple wrapper for protein sequence data."""
    def __init__(self, data: List[Tuple[str, str]]):
        self.data = data 

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]

class ESMFactory:
    """A memory-efficient pipeline for embedding thousands of protein sequences using ESM-2."""
    
    ESM_MODELS_DICT = {
        320: "esm2_t6_8M_UR50D",
        480: "esm2_t12_35M_UR50D",
        640: "esm2_t30_150M_UR50D",
        1280: "esm2_t33_650M_UR50D",
        2560: "esm2_t36_3B_UR50D",
        5120: "esm2_t48_15B_UR50D"
    }

    def __init__(self, embedding_size: int = 1280, device: str = None):
        if embedding_size not in self.ESM_MODELS_DICT:
            raise ValueError(f"Invalid size. Choose from: {list(self.ESM_MODELS_DICT.keys())}")

        self.device = device or ('cuda' if torch.cuda.is_available() else 'cpu')
        model_name = self.ESM_MODELS_DICT[embedding_size]
        
        print(f"Loading {model_name} to {self.device}...")
        self.model, self.alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        self.model.to(self.device).eval()
        
        self.batch_converter = self.alphabet.get_batch_converter()
        self.repr_layer = self.model.num_layers
        self.max_tokens = 1022 

    def _get_dataloader(self, sequence_data: List[Tuple[str, str]], batch_size: int) -> DataLoader:
        """Creates a DataLoader with the ESM batch converter."""
        dataset = ProteinDataset(sequence_data)
        return DataLoader(dataset, batch_size=batch_size, shuffle=False, collate_fn=self.batch_converter)

    def _process_batches(self, dataloader: DataLoader, desc: str) -> Iterator[Tuple[str, torch.Tensor]]:
        """
        Core generator that handles the forward pass, truncation, and padding removal.
        Yields (sequence_id, residue_level_tensor_on_cpu) to prevent memory bloat.
        """
        with torch.no_grad():
            for labels, _, tokens in dataloader:
                # Truncate to max model length
                tokens = tokens[:, :self.max_tokens].to(self.device)
                
                # Forward pass
                out = self.model(tokens, repr_layers=[self.repr_layer])
                token_embeddings = out["representations"][self.repr_layer]
                
                # Calculate true lengths
                batch_lens = (tokens != self.alphabet.padding_idx).sum(1)
                
                for i, length in enumerate(batch_lens):
                    seq_id = labels[i]
                    # Slice out BOS (0) and EOS (length-1), move to CPU
                    residue_reps = token_embeddings[i, 1 : int(length) - 1].cpu()
                    
                    yield seq_id, residue_reps

    def get_mean_embeddings(self, sequence_data: List[Tuple[str, str]], batch_size: int = 32) -> Dict[str, torch.Tensor]:
        """Returns a dictionary of {id: mean_vector_tensor} (Shape: D)"""
        dataloader = self._get_dataloader(sequence_data, batch_size)
        
        # Uses dictionary comprehension to calculate the mean on the fly
        return {
            seq_id: reps.mean(dim=0)
            for seq_id, reps in self._process_batches(dataloader, "Computing Mean Embeddings")
        }

    def get_per_residue_embeddings(self, sequence_data: List[Tuple[str, str]], batch_size: int = 16) -> Dict[str, torch.Tensor]:
        """Returns a dictionary of {id: residue_tensor} (Shape: L x D)"""
        dataloader = self._get_dataloader(sequence_data, batch_size)
        
        # Directly constructs the dictionary from the generator
        return {
            seq_id: reps 
            for seq_id, reps in self._process_batches(dataloader, "Computing Residue Embeddings")
        }

    def get_single_sequence_mean(self, sequence: str) -> torch.Tensor:
        """Helper for a single string input."""
        return self.get_mean_embeddings([("query", sequence)], batch_size=1)["query"]

# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(path: Path) -> List[Tuple[str, str]]:
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


# ---------------------------------------------------------------------------
# Checkpointed embedding
# ---------------------------------------------------------------------------

def embed_with_checkpoint(factory: ESMFactory,
                           records: List[Tuple[str, str]],
                           output_path: Path,
                           batch_size: int) -> Dict[str, torch.Tensor]:
    """
    Run get_mean_embeddings with per-batch checkpointing.
    Resumes from checkpoint if one exists.
    Saves final .pt and removes checkpoint on success.
    """
    checkpoint_path = output_path.with_suffix(".checkpoint.pt")

    # Resume from checkpoint if available
    done: Dict[str, torch.Tensor] = {}
    if checkpoint_path.exists():
        done = torch.load(checkpoint_path, map_location="cpu")
        log.info(f"  Resuming: {len(done):,} sequences already embedded")

    todo = [(h, s) for h, s in records if h not in done]
    log.info(f"  To embed: {len(todo):,} sequences")

    if todo:
        # Batch manually so we can checkpoint after each batch
        for start in range(0, len(todo), batch_size):
            batch = todo[start:start + batch_size]
            batch_embeddings = factory.get_mean_embeddings(batch, batch_size=len(batch))
            done.update(batch_embeddings)
            torch.save(done, checkpoint_path)
            log.info(f"  Checkpoint: {len(done):,}/{len(records):,} sequences done")

    # Save final output
    torch.save(done, output_path)
    log.info(f"  Saved embeddings -> {output_path}")

    if checkpoint_path.exists():
        checkpoint_path.unlink()
        log.info("  Checkpoint removed")

    return done


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute ESM-2 mean embeddings for a FASTA file."
    )
    parser.add_argument(
        "--fasta", type=Path, required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output .pt file path"
    )
    parser.add_argument(
        "--embedding-size", type=int, default=1280,
        choices=[320, 480, 640, 1280, 2560, 5120],
        help="ESM-2 embedding dimension (selects model, default: 1280 = 650M params)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=32,
        help="Sequences per batch (reduce if OOM, default: 32)"
    )
    parser.add_argument(
        "--device", type=str, default=None,
        help="torch device (default: auto-detect cuda/cpu)"
    )
    args = parser.parse_args()

    if not args.fasta.exists():
        log.error(f"FASTA not found: {args.fasta}")
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.output.exists():
        log.info(f"Output already exists, skipping: {args.output}")
        sys.exit(0)

    # Log device info
    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"Device: {device}")
    if device == "cuda":
        log.info(f"  GPU: {torch.cuda.get_device_name(0)}")
        log.info(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    records = parse_fasta(args.fasta)
    factory = ESMFactory(embedding_size=args.embedding_size, device=args.device)

    embeddings = embed_with_checkpoint(factory, records, args.output, args.batch_size)

    log.info(f"\nDone. {len(embeddings):,} embeddings saved to {args.output}")
    log.info(f"Embedding shape: {next(iter(embeddings.values())).shape}")


if __name__ == "__main__":
    main()