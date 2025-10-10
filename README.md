# Tools for the Sticks and Lines Project

This project is dockerized and employs a devcontainer with VSCode for reproducibility. If Docker is installed, the devcontainer can be spun up by following the steps below.

## Prerequisites

- [Docker Desktop](https://www.docker.com/products/docker-desktop) installed and running
- [Visual Studio Code](https://code.visualstudio.com/) installed
- [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) for VSCode

## Getting Started

1. **Open the project in VSCode**
   ```bash
   code /path/to/sticks-and-lines-project
   ```

2. **Reopen in Container**
   - When prompted, click "Reopen in Container"
   - Or use Command Palette (Ctrl+Shift+P / Cmd+Shift+P) and select `Dev Containers: Reopen in Container`


## Project Structure

```
├── README.md                      # This file
├── checksum.txt                   # Data integrity checksums
├── data/                          # Data directory
│   ├── ObjectOrientation_3_0.csv
│   ├── PrefViews All 3D objects
│   ├── PrefViews Null 3D objects
│   ├── experimental_datasets/     # Experimental data
│   │   └── PVMDA 7-1 Dataset/
│   │       ├── PVMDA_combined_rotation_and_moments.csv
│   │       ├── annotated_images.zip
│   │       ├── intermediate_csvs/
│   │       ├── renditions.zip
│   │       └── second_order_histograms.zip
│   └── null_datasets/             # Null hypothesis datasets
│       └── MT19937/
│           ├── MT19937_combined_rotation_and_moments.csv
│           ├── annotated_images.zip
│           ├── intermediate_csvs/
│           ├── renditions.zip
│           ├── runtime.log
│           └── second_order_histograms.zip
├── lib/                           # External libraries
│   └── image_entropy/             # MATLAB entropy calculation functions
│       ├── cpda.m
│       ├── edge_entropy.m
│       ├── scalarq.m
│       ├── shannon_entropy.m
│       └── stheta.m
├── requirements.txt               # Python dependencies
└── src/                           # Source code
    ├── ellipse_orientation_script.m
    ├── fit_ellipse_conv_hull.m
    ├── fit_ellipse_filled_region.m
    ├── fit_ellipse_pure_moments_cov.m
    ├── generate_experimental_renditions.ipynb
    ├── generate_null_rotations.ipynb
    ├── investigate_image_discrepancy.ipynb
    ├── major_axis_features.py
    ├── region_props_multi.m
    ├── simple_region_props.m
    └── test_main_axis_features.ipynb
```

## Key Components

### Data Processing
- **experimental_datasets/**: Contains PVMDA experimental data with rotation and moment calculations
- **null_datasets/**: MT19937-generated null hypothesis data for comparison
- Intermediate CSVs include region properties, histograms, and entropy measurements

### Analysis Scripts
- **MATLAB scripts** (`.m`): Ellipse fitting and region property calculations
- **Python notebooks** (`.ipynb`): Data generation and analysis workflows
- **major_axis_features.py**: Python module for feature extraction

### Libraries
- **image_entropy/**: External MATLAB functions for various entropy calculations (Shannon, edge, etc.)

## Usage

After the devcontainer is running, you can:

- Run Jupyter notebooks in the `src/` directory
- Execute MATLAB scripts for region property analysis
- Process experimental and null datasets
- Generate renditions and analyze object orientations

## Troubleshooting

- **Container won't start**: Ensure Docker Desktop is running
- **Missing dependencies**: Rebuild the container using `Dev Containers: Rebuild Container`

## Notes

- All paths are relative to the project root
- Generated files and caches are excluded from version control
- Runtime logs are saved in dataset directories for reproducibility