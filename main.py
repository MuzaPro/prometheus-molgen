#!/usr/bin/env python3
"""
Prometheus MolGen: Organic Chemistry Molecule SVG Generator
Main entry point for the application.
"""

import os
import sys
import click
from typing import List, Optional

# Import core modules
from src.molecule_parser import parse_molecule
from src.svg_generator import generate_svg
from src.file_manager import save_svg, read_molecule_list

# Define version
__version__ = "0.1.0"

@click.group()
@click.version_option(__version__)
def cli():
    """Prometheus MolGen: Generate SVG representations of organic chemistry molecules."""
    pass

@cli.command()
@click.option('--smiles', help='SMILES notation of a molecule')
@click.option('--name', help='Name of the molecule (used for filename)')
@click.option('--input', help='Input file with molecule definitions')
@click.option('--output', default='output', help='Output directory for SVG files')
@click.option('--style', default='modern', help='Styling to apply (modern, classic)')
def generate(smiles: Optional[str], name: Optional[str], input: Optional[str], 
             output: str, style: str):
    """Generate SVG representations of molecules."""
    
    # Create output directory if it doesn't exist
    os.makedirs(output, exist_ok=True)
    
    # Process a single molecule if provided via SMILES
    if smiles and name:
        click.echo(f"Generating SVG for {name} from SMILES: {smiles}")
        molecule = parse_molecule(smiles)
        if molecule:
            svg_content = generate_svg(molecule, style=style)
            save_svg(svg_content, name, output_dir=output)
            click.echo(f"Generated: {output}/{name}.svg")
        else:
            click.echo(f"Error: Could not parse molecule: {smiles}", err=True)
        return
    
    # Process molecules from input file
    if input:
        click.echo(f"Processing molecules from: {input}")
        try:
            molecules = read_molecule_list(input)
            
            with click.progressbar(molecules, label='Generating SVGs') as bar:
                for mol_name, mol_smiles in bar:
                    molecule = parse_molecule(mol_smiles)
                    if molecule:
                        svg_content = generate_svg(molecule, style=style)
                        save_svg(svg_content, mol_name, output_dir=output)
                    else:
                        click.echo(f"Warning: Could not parse molecule: {mol_name}", err=True)
            
            click.echo(f"Successfully generated SVGs in: {output}")
        except Exception as e:
            click.echo(f"Error processing input file: {e}", err=True)
            return
    
    # If no input specified
    if not smiles and not input:
        click.echo("Error: Please provide either --smiles and --name or --input", err=True)

@cli.command()
@click.option('--input', help='SVG file to validate', required=True)
def validate(input: str):
    """Validate SVG files against specifications."""
    # This is a placeholder for future validation functionality
    click.echo(f"Validation for {input} not yet implemented")

if __name__ == '__main__':
    cli()
