"""
SVG Generator Module
Creates SVG representations of molecules with custom styling
"""

import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D
from typing import Dict, Any, Optional, List, Tuple
import io
import re
import math
import svgwrite

from src.molecule_parser import get_molecule_data
from src.style_manager import get_style_settings

def generate_svg(molecule: Chem.Mol, style: str = 'modern', 
                 size: Tuple[int, int] = (100, 100)) -> str:
    """
    Generate an SVG representation of a molecule with custom styling.
    
    Args:
        molecule: RDKit molecule object
        style: Styling to apply ('modern', 'classic')
        size: Size of the SVG output (width, height)
        
    Returns:
        SVG content as a string
    """
    # First get the molecule data
    mol_data = get_molecule_data(molecule)
    
    # Apply the selected style
    style_settings = get_style_settings(style)
    
    # Generate initial SVG using RDKit
    initial_svg = generate_initial_svg(molecule, size)
    
    # Apply custom styling to the SVG
    styled_svg = apply_custom_styling(initial_svg, mol_data, style_settings)
    
    return styled_svg

def generate_initial_svg(molecule: Chem.Mol, size: Tuple[int, int] = (100, 100)) -> str:
    """
    Generate initial SVG using RDKit's drawing capabilities.
    
    Args:
        molecule: RDKit molecule object
        size: Size of the SVG output (width, height)
        
    Returns:
        Initial SVG content as a string
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    
    # Clean up SVG (remove RDKit specific elements if needed)
    svg = svg.replace('xmlns:rdkit="http://www.rdkit.org/xml"', '')
    svg = svg.replace('xmlns:xlink="http://www.w3.org/1999/xlink"', '')
    
    return svg

def apply_custom_styling(svg: str, mol_data: Dict[str, Any], 
                         style_settings: Dict[str, Any]) -> str:
    """
    Apply custom styling to the SVG based on style settings.
    
    Args:
        svg: Initial SVG content
        mol_data: Molecule data from get_molecule_data()
        style_settings: Style settings from get_style_settings()
        
    Returns:
        Styled SVG content as a string
    """
    # Parse the SVG
    try:
        svg_tree = ET.fromstring(svg)
    except ET.ParseError:
        # Add missing root element if necessary
        svg = f"<svg>{svg}</svg>"
        svg_tree = ET.fromstring(svg)
    
    # Set viewBox and dimensions
    svg_tree.set('width', str(style_settings['svg_width']))
    svg_tree.set('height', str(style_settings['svg_height']))
    svg_tree.set('viewBox', f"0 0 {style_settings['svg_width']} {style_settings['svg_height']}")
    
    # Find and style all atom text elements
    for text_elem in svg_tree.findall(".//text"):
        # Apply font styling
        text_elem.set('font-family', style_settings['font_family'])
        text_elem.set('font-size', str(style_settings['atom_font_size']))
        text_elem.set('font-weight', style_settings['font_weight'])
        text_elem.set('text-anchor', 'middle')
        text_elem.set('dominant-baseline', 'central')
        
        # Style by atom type if needed
        atom_symbol = text_elem.text.strip() if text_elem.text else ""
        if atom_symbol in style_settings['atom_colors']:
            text_elem.set('fill', style_settings['atom_colors'][atom_symbol])
        else:
            text_elem.set('fill', style_settings['default_atom_color'])
    
    # Style bonds
    for path_elem in svg_tree.findall(".//path"):
        # Check if this is a bond path
        if 'class' in path_elem.attrib and 'bond' in path_elem.attrib['class']:
            path_elem.set('stroke', style_settings['bond_color'])
            path_elem.set('stroke-width', str(style_settings['bond_stroke_width']))
    
    # Convert back to string
    return ET.tostring(svg_tree, encoding='unicode')

def create_svg_from_scratch(mol_data: Dict[str, Any], 
                           style_settings: Dict[str, Any]) -> str:
    """
    Create an SVG from scratch based on molecule data and style settings.
    This is an alternative to modifying RDKit's output and gives more control.
    
    Args:
        mol_data: Molecule data from get_molecule_data()
        style_settings: Style settings from get_style_settings()
        
    Returns:
        SVG content as a string
    """
    # Calculate scaling to fit molecule in viewBox
    width, height = mol_data['dimensions']['width'], mol_data['dimensions']['height']
    min_x, min_y = mol_data['dimensions']['min_x'], mol_data['dimensions']['min_y']
    
    # Create aspect-ratio preserving scaling
    max_dim = max(width, height)
    scale = (style_settings['svg_width'] * 0.8) / max_dim  # Use 80% of available space
    
    # Create SVG drawing
    dwg = svgwrite.Drawing(profile='tiny', 
                          size=(style_settings['svg_width'], style_settings['svg_height']))
    
    # Set viewBox
    dwg.viewbox(0, 0, style_settings['svg_width'], style_settings['svg_height'])
    
    # Create a group for the molecule, centered in the viewBox
    mol_group = dwg.g(transform=f"translate({style_settings['svg_width']/2 - width*scale/2}, "
                               f"{style_settings['svg_height']/2 - height*scale/2})")
    
    # Draw bonds
    for bond in mol_data['bonds']:
        begin_atom = next(a for a in mol_data['atoms'] if a['idx'] == bond['begin_atom_idx'])
        end_atom = next(a for a in mol_data['atoms'] if a['idx'] == bond['end_atom_idx'])
        
        # Calculate scaled coordinates
        x1, y1 = (begin_atom['x'] - min_x) * scale, (begin_atom['y'] - min_y) * scale
        x2, y2 = (end_atom['x'] - min_x) * scale, (end_atom['y'] - min_y) * scale
        
        # Determine bond type and draw accordingly
        bond_type = str(bond['bond_type'])
        
        if bond_type == 'SINGLE':
            line = dwg.line(start=(x1, y1), end=(x2, y2), 
                          stroke=style_settings['bond_color'], 
                          stroke_width=style_settings['bond_stroke_width'])
            mol_group.add(line)
        
        elif bond_type == 'DOUBLE':
            # Calculate parallel line offset
            angle = math.atan2(y2 - y1, x2 - x1) + math.pi/2
            offset = style_settings['double_bond_spacing'] / 2
            dx, dy = math.cos(angle) * offset, math.sin(angle) * offset
            
            line1 = dwg.line(start=(x1 + dx, y1 + dy), end=(x2 + dx, y2 + dy), 
                           stroke=style_settings['bond_color'], 
                           stroke_width=style_settings['bond_stroke_width'])
            line2 = dwg.line(start=(x1 - dx, y1 - dy), end=(x2 - dx, y2 - dy), 
                           stroke=style_settings['bond_color'], 
                           stroke_width=style_settings['bond_stroke_width'])
            
            mol_group.add(line1)
            mol_group.add(line2)
        
        elif bond_type == 'TRIPLE':
            # Center line
            line1 = dwg.line(start=(x1, y1), end=(x2, y2), 
                           stroke=style_settings['bond_color'], 
                           stroke_width=style_settings['bond_stroke_width'])
            
            # Calculate parallel line offset
            angle = math.atan2(y2 - y1, x2 - x1) + math.pi/2
            offset = style_settings['triple_bond_spacing']
            dx, dy = math.cos(angle) * offset, math.sin(angle) * offset
            
            line2 = dwg.line(start=(x1 + dx, y1 + dy), end=(x2 + dx, y2 + dy), 
                           stroke=style_settings['bond_color'], 
                           stroke_width=style_settings['bond_stroke_width'])
            line3 = dwg.line(start=(x1 - dx, y1 - dy), end=(x2 - dx, y2 - dy), 
                           stroke=style_settings['bond_color'], 
                           stroke_width=style_settings['bond_stroke_width'])
            
            mol_group.add(line1)
            mol_group.add(line2)
            mol_group.add(line3)
    
    # Draw atoms
    for atom in mol_data['atoms']:
        x, y = (atom['x'] - min_x) * scale, (atom['y'] - min_y) * scale
        
        # Get atom color
        atom_color = style_settings['atom_colors'].get(
            atom['symbol'], style_settings['default_atom_color'])
        
        # Create atom text with proper styling
        atom_text = dwg.text(atom['symbol'], 
                           insert=(x, y),
                           font_family=style_settings['font_family'],
                           font_size=style_settings['atom_font_size'],
                           font_weight=style_settings['font_weight'],
                           fill=atom_color,
                           text_anchor='middle',
                           dominant_baseline='central')
        
        mol_group.add(atom_text)
        
        # Add charge if not zero
        if atom['charge'] != 0:
            charge_text = '+' if atom['charge'] > 0 else '-'
            if abs(atom['charge']) > 1:
                charge_text = f"{abs(atom['charge'])}{charge_text}"
                
            charge = dwg.text(charge_text,
                            insert=(x + style_settings['charge_offset_x'], 
                                    y + style_settings['charge_offset_y']),
                            font_family=style_settings['font_family'],
                            font_size=style_settings['charge_font_size'],
                            fill=atom_color,
                            text_anchor='middle')
            
            mol_group.add(charge)
    
    # Add the molecule group to the drawing
    dwg.add(mol_group)
    
    # Return the SVG as a string
    return dwg.tostring()
