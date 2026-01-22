# 3D-GRID-SPLITTER

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

**Automated Grid Slicing Tool for Large 3D Models**

This repository provides a Python-based solution to split large 3D mesh files (`.stl`, `.obj`) into smaller, user-defined grid chunks (e.g., 50mm x 50mm cubes).

Unlike simple cutting tools, this script automatically **fills the cut surfaces (capping)**, making the resulting parts "watertight" and immediately ready for **3D printing** or **segment-based analysis**.

## 🚀 Features

- **Grid Slicing:** Automatically splits a 3D model into an `N x M` grid based on specified dimensions (mm).
- **Watertight Capping:** Automatically fills the cut surfaces (closes holes) to ensure the parts are solid.
- **Visualization:** Includes a built-in 3D viewer (PyVista) to preview the split chunks.
- **Format Support:** Supports `.stl`, `.obj`, and other formats supported by PyVista.

## 📦 Installation

1. Clone the repository:
   ```bash
   git clone [https://github.com/YEONCHEOL-HA/3D-GRID-SPLITTER.git](https://github.com/YEONCHEOL-HA/3D-GRID-SPLITTER.git)
   cd 3D-GRID-SPLITTER

<img width="1024" height="768" alt="image" src="https://github.com/user-attachments/assets/d9ce6541-015e-4e8c-ac53-84655c0406fd" />
