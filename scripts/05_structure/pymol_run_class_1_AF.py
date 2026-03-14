from pymol import cmd
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

STRUCTURE_DIR = Path(PATHS["results"]) / "structure"
STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)


def setup_fuma_scene():
    cmd.reinitialize()

    # 1. Load AlphaFold model for FumA (P0AC33)
    cmd.load("https://alphafold.ebi.ac.uk/files/AF-P0AC33-F1-model_v6.pdb", "fuma")

    # 2. Visual styling
    cmd.bg_color("white")
    cmd.color("palegreen", "fuma")
    cmd.set("cartoon_fancy_helices", 1)

    # 3. Active site: Cys coordination for Fe-S cluster
    cmd.select("active_site", "fuma and resi 105+224+318")
    cmd.show("spheres", "active_site")
    cmd.set("sphere_scale", 0.6, "active_site")
    cmd.color("brightorange", "active_site")
    cmd.util.cnc("active_site")

    # 4. Surface Highlight for the Pocket
    cmd.show("surface", "active_site")
    cmd.set("transparency", 0.4, "active_site")

    # 5. Publication-Level Rendering Settings
    cmd.set("ambient", 0.5)
    cmd.set("ray_shadows", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_opaque_background", "on")

    # 6. Adjust view
    cmd.orient("fuma")
    cmd.zoom("all", buffer=-15)

    # 7. Export image
    export_path = STRUCTURE_DIR / "FumA_ClassI.png"
    print("Rendering high-resolution image...")
    cmd.ray(1200, 800)
    cmd.png(str(export_path), dpi=600)
    print(f"Image saved to: {export_path}")

    # 8. Save session
    session_path = STRUCTURE_DIR / "FumA_ClassI_analysis.pse"
    cmd.save(str(session_path))
    print(f"PyMOL session saved to: {session_path}")
    print("FumA render complete — Cys105, Cys224, Cys318.")


setup_fuma_scene()