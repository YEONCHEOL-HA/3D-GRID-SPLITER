import pyvista as pv
import numpy as np
import os
import argparse
import sys

def split_mesh(input_path, output_dir, grid_size, close_hole=True, visualize=False):
    """
    Splits a 3D mesh into a grid of smaller chunks.
    """
    
    # 1. Load Mesh
    if not os.path.exists(input_path):
        print(f"[Error] File not found: {input_path}")
        return

    print(f"Loading mesh: {input_path}...")
    mesh = pv.read(input_path)
    
    # 2. Create Output Directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # 3. Get Bounds and Set Ranges
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
    gx, gy, gz = grid_size
    
    print(f"Original Bounds: X({xmin:.1f}~{xmax:.1f}), Y({ymin:.1f}~{ymax:.1f}), Z({zmin:.1f}~{zmax:.1f})")
    print(f"Grid Size: {gx}mm x {gy}mm x {gz}mm")

    x_ranges = np.arange(xmin, xmax, gx)
    y_ranges = np.arange(ymin, ymax, gy)
    z_ranges = np.arange(zmin, zmax, gz)

    # 4. Visualization Setup
    p = None
    if visualize:
        # Check for headless environment (e.g., Colab, Server)
        if 'google.colab' in sys.modules:
            pv.start_xvfb()
        
        p = pv.Plotter()
        p.add_mesh(mesh, opacity=0.1, color='lightgrey', style='wireframe')

    count = 0
    total_steps = len(x_ranges) * len(y_ranges) * len(z_ranges)
    print(f"Processing split... (Max chunks: {total_steps})")
    
    # 5. Split Loop
    for i, x_start in enumerate(x_ranges):
        for j, y_start in enumerate(y_ranges):
            for k, z_start in enumerate(z_ranges):
                
                # Define clipping bounds
                x_end = min(x_start + gx, xmax)
                y_end = min(y_start + gy, ymax)
                z_end = min(z_start + gz, zmax)
                clip_bounds = [x_start, x_end, y_start, y_end, z_start, z_end]
                
                try:
                    # Clip the mesh
                    clipped = mesh.clip_box(clip_bounds, invert=False)
                    
                    if clipped.n_cells > 0:
                        # Extract surface
                        surface = clipped.extract_surface()
                        
                        # Close holes (Cap the cut surfaces)
                        if close_hole:
                            surface = surface.fill_holes(hole_size=10000000)
                            surface = surface.triangulate()

                        # Save file
                        filename = f"part_x{i}_y{j}_z{k}.stl"
                        filepath = os.path.join(output_dir, filename)
                        surface.save(filepath)
                        count += 1
                        print(f" - Saved: {filename} (Cells: {surface.n_cells})")

                        # Add to visualization
                        if visualize and p:
                            p.add_mesh(surface, color=np.random.random(3))
                            p.add_mesh(pv.Box(clip_bounds), style='wireframe', color='black', opacity=0.2)
                            
                except Exception as e:
                    print(f"[Warning] Error at index {i},{j},{k}: {e}")
                    continue

    print(f"\nDone! Total {count} parts saved in '{output_dir}'.")
    
    if visualize and p and count > 0:
        print("Opening visualization window...")
        p.show_grid()
        p.show_axes()
        # Use 'static' for compatibility, or remove arguments for interactive window
        p.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a 3D Mesh (STL/OBJ) into smaller grid chunks.")
    
    parser.add_argument("--input", "-i", type=str, required=True, help="Path to input 3D file (.stl, .obj)")
    parser.add_argument("--output", "-o", type=str, default="output_parts", help="Directory to save split files")
    parser.add_argument("--size", "-s", type=float, nargs=3, default=[50, 50, 50], help="Grid size in mm (x y z). Default: 50 50 50")
    parser.add_argument("--no-cap", action="store_true", help="Disable hole filling (capping)")
    parser.add_argument("--vis", action="store_true", help="Enable visualization of the result")

    args = parser.parse_args()
    
    # Invert the no-cap logic for the function
    close_holes = not args.no_cap
    
    split_mesh(args.input, args.output, args.size, close_holes, args.vis)
