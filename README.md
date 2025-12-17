# GLI Structure Analysis Project

This project contains tools and scripts for analyzing GLI protein structures, including molecular docking analysis with GANT61 and structural validation tools.

## Environment Setup

### Creating Anaconda Virtual Environment

To set up the required environment for this project, create a new Anaconda virtual environment with the necessary dependencies:

```bash
# Create a new conda environment named 'gli-docking' with Python 3.9
conda create -n gli-docking python=3.9

# Activate the environment
conda activate gli-docking

# Install required packages
conda install -c conda-forge biopython numpy pandas matplotlib scipy
conda install -c conda-forge rdkit pymol-open-source
pip install prody MDAnalysis

# For molecular docking (optional)
conda install -c conda-forge openbabel
```

### Installing Missing Packages in Existing Environment

If you already have a conda environment (like `gli-docking`) and need to add missing packages:

```bash
# Activate your existing environment
conda activate gli-docking

# Install rdkit and other missing packages
conda install -c conda-forge rdkit
conda install -c conda-forge biopython numpy pandas matplotlib scipy
conda install -c conda-forge pymol-open-source openbabel
pip install prody MDAnalysis
```

### Alternative: Using environment.yml

You can also create the environment using the provided environment file:

```bash
conda env create -f environment.yml
conda activate gli-docking
```

## Project Structure

- **PDB Files**: Contains GLI protein structures and GANT61 ligand files
- **Analysis Scripts**: Python scripts for structural analysis and binding site validation
- **Docking Tools**: Scripts for molecular docking preparation and analysis

## Usage

After activating the conda environment, you can run the analysis scripts:

```bash
python analyze_interactions.py
python binding_site_validator.py
python gant61_binding_analysis.py
```

Make sure to activate the conda environment before running any scripts to ensure all dependencies are available.
