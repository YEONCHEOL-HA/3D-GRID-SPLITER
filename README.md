# 3D Mesh Grid Splitter

A Python tool to split large 3D meshes (STL, OBJ) into smaller, uniform grid chunks. This is particularly useful for **3D printing oversized models** on smaller printers or for **spatial data analysis**.

## Features

- **Grid Slicing:** Automatically splits a 3D model into an `N x M` grid based on specified dimensions (mm).
- **Watertight Capping:** Automatically fills the cut surfaces (closes holes) to ensure the parts are solid and 3D printable.
- **Visualization:** Provides a 3D preview of the split parts with a wireframe grid.
- **Format Support:** Supports `.stl`, `.obj`, and other formats supported by PyVista.

## Requirements

- Python 3.7+
- `pyvista`
- `numpy`

## Installation

1. Clone the repository:
   ```bash
   git clone [https://github.com/YOUR_USERNAME/3d-mesh-grid-splitter.git](https://github.com/YOUR_USERNAME/3d-mesh-grid-splitter.git)
   cd 3d-mesh-grid-splitter
