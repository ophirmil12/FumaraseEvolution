from pymol import cmd
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

STRUCTURE_DIR = Path(PATHS["results"]) / "structure"
STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)


def setup_fumc_scene():
    cmd.reinitialize()

    # 1. Load AlphaFold model for FumC (P05042)
    cmd.load("https://alphafold.ebi.ac.uk/files/AF-P05042-F1-model_v6.pdb", "fumc")

    # 2. Visual styling
    cmd.bg_color("white")
    cmd.color("palegreen", "fumc")
    cmd.set("cartoon_fancy_helices", 1)

    # 3. Active site
    cmd.select("active_site", "fumc and resi 188+318+324+331")
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
    cmd.orient("fumc")
    cmd.zoom("all", buffer=-5)

    # 7. Export image
    export_path = STRUCTURE_DIR / "FumC_ClassII.png"
    print("Rendering high-resolution image...")
    cmd.ray(1200, 800)
    cmd.png(str(export_path), dpi=600)
    print(f"Image saved to: {export_path}")

    # 8. Save session
    session_path = STRUCTURE_DIR / "FumC_ClassII_analysis.pse"
    cmd.save(str(session_path))
    print(f"PyMOL session saved to: {session_path}")
    print("FumC render complete — His188, Ser318, Lys324, Glu331.")


setup_fumc_scene()