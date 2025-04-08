# Getting Started with Prometheus MolGen

This guide will help you set up and start using Prometheus MolGen to generate SVG representations of organic chemistry molecules.

## Installation

### Prerequisites

Before installing Prometheus MolGen, make sure you have:

1. **Python 3.7 or higher** installed on your system
2. **Git** installed (for cloning the repository)
3. Basic familiarity with command-line operations

### Step 1: Clone the Repository

```bash
git clone https://github.com/your-username/prometheus-molgen.git
cd prometheus-molgen
```

### Step 2: Set Up a Virtual Environment

It's recommended to use a virtual environment to avoid package conflicts:

```bash
# Create a virtual environment
python -m venv venv

# Activate the virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate
```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

> **Note**: RDKit can sometimes be challenging to install via pip. If you encounter issues, consider using conda: `conda install -c conda-forge rdkit`

## Basic Usage

### Generating a Single Molecule SVG

To generate an SVG for a single molecule using SMILES notation:

```bash
python main.py generate --smiles "CO" --name "methanol"
```

This will create `methanol.svg` in the default `output` directory.

### Processing Multiple Molecules

To generate SVGs for multiple molecules defined in a file:

```bash
python main.py generate --input examples/molecules.csv --output assets/molecules/
```

This will process all molecules in `molecules.csv` and save the SVGs to the `assets/molecules/` directory.

## Input Format Examples

### CSV Format

Create a CSV file with the following structure:

```csv
name,smiles
Methanol,CO
Ethanol,CCO
Benzene,c1ccccc1
```

### Text Format

Alternatively, create a plain text file with space-separated values:

```
Methanol CO
Ethanol CCO
Benzene c1ccccc1
```

### JSON Format

Or use JSON for more flexibility:

```json
[
  {"name": "Methanol", "smiles": "CO"},
  {"name": "Ethanol", "smiles": "CCO"},
  {"name": "Benzene", "smiles": "c1ccccc1"}
]
```

## Styling Options

To apply a specific style to your generated SVGs:

```bash
python main.py generate --smiles "c1ccccc1" --name "benzene" --style modern
```

Available styles:
- `modern` (default): Larger atom labels, shorter bonds
- `classic`: Traditional chemical drawing style
- `compact`: Space-efficient style for complex molecules

## Integration with the Quiz Game

To generate SVGs directly for your Organic Chemistry Quiz Game:

1. Make sure both repositories are cloned on your system
2. Run the following command:

```bash
python main.py generate --input molecules.csv --output ../Organic-Chemistry-Quiz-JS/assets/molecules/
```

3. Update the paths in your quiz game's `questions.json` file to reference these SVGs:

```json
"product": {
  "name": "Methyl Ethyl Ether",
  "formula": "CH3OCH2CH3",
  "imagePath": "assets/molecules/methyl-ethyl-ether.svg",
  "textRepresentation": "CH3-O-CH2CH3"
}
```

## Common SMILES Notations

Here are some common SMILES notations for reference:

| Molecule | SMILES |
|----------|--------|
| Methanol | CO |
| Ethanol | CCO |
| Acetic Acid | CC(=O)O |
| Benzene | c1ccccc1 |
| Methyl Ethyl Ether | CCOC |
| tert-Butyl Chloride | CC(C)(C)Cl |
| Sodium Hydroxide | [Na+].[OH-] |
| Water | O |

## Troubleshooting

### RDKit Installation Issues

If you're having trouble installing RDKit:

1. Try using conda: `conda install -c conda-forge rdkit`
2. Check the [RDKit documentation](https://www.rdkit.org/docs/Install.html) for your platform
3. Consider using a pre-built Docker image with RDKit installed

### SVG Generation Problems

If molecules aren't displaying correctly:

1. Check that your SMILES notation is correct
2. Try generating with the `classic` style to see if it's a styling issue
3. For complex molecules, try the `compact` style

## Next Steps

- Experiment with different molecules and styles
- Consider contributing to the project if you identify improvements
- Check out the project's GitHub issues for feature requests or bugs you might help with

For more detailed information, refer to the full [README.md](../README.md) file.