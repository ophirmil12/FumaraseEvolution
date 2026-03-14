from pymol import cmd
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

STRUCTURE_DIR = Path(PATHS["results"]) / "structure"
STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)


def setup_1YFE_scene():
    cmd.reinitialize()

    # 1. Fetch and setup states
    cmd.fetch("1YFE")
    cmd.set("all_states", "on")
    cmd.split_states("1YFE")

    # 2. Generate symmetry neighbors
    cmd.symexp("sym", "1YFE", "(1YFE)", 4)

    # 3. Cleanup: Remove solvent
    cmd.remove("solvent")

    # 4. Filter objects
    keep_list = ["1YFE", "sym01000000", "sym02000000", "sym03000000"]
    all_objs = cmd.get_names("objects")
    for obj in all_objs:
        if obj not in keep_list:
            cmd.delete(obj)

    # 5. Visual styling: White background and Pastel colors
    cmd.bg_color("white")
    colors = ["lightmagenta", "lightblue", "palegreen", "wheat"]
    for i, obj in enumerate(keep_list):
        if i < len(colors):
            cmd.color(colors[i], obj)

    # 6. Active site: High-visibility styling
    cmd.select("active_site", "resi 187+317+323+330")
    cmd.show("spheres", "active_site")
    cmd.set("sphere_scale", 0.6, "active_site")
    cmd.color("brightorange", "active_site")
    cmd.util.cnc("active_site")

    # 7. Surface Highlight for the Pocket
    cmd.show("surface", "active_site")
    cmd.set("transparency", 0.4, "active_site")

    # 8. Publication-Level Rendering Settings
    cmd.set("ambient", 0.5)
    cmd.set("ray_shadows", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_opaque_background", "on")

    cmd.zoom("all", buffer=2)

    # 9. Save session
    session_path = STRUCTURE_DIR / "1YFE_analysis.pse"
    cmd.save(str(session_path))
    print(f"PyMOL session saved to: {session_path}")

    # 10. Export image
    export_path = STRUCTURE_DIR / "1YFE_final_render.png"
    print("Rendering high-resolution image...")
    cmd.ray(1200, 800)
    cmd.png(str(export_path), dpi=600)
    print(f"Image saved to: {export_path}")


setup_1YFE_scene()