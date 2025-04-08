"""
Style Manager Module
Manages styling rules and templates for SVG generation
"""

from typing import Dict, Any

def get_style_settings(style_name: str = 'modern') -> Dict[str, Any]:
    """
    Get the style settings for the specified style.
    
    Args:
        style_name: Name of the style ('modern', 'classic')
        
    Returns:
        Dictionary of style settings
    """
    # Base settings common to all styles
    base_settings = {
        'svg_width': 100,
        'svg_height': 100,
        'default_atom_color': 'black',
        'bond_color': 'black',
        'font_family': 'Arial, sans-serif',
        'font_weight': 'normal',
        'atom_colors': {
            'C': 'black',
            'H': 'black',
            'O': 'black',
            'N': 'black',
            'F': 'black',
            'Cl': 'black',
            'Br': 'black',
            'I': 'black',
            'P': 'black',
            'S': 'black',
            'Na': 'black',
            'Mg': 'black',
            'K': 'black'
        },
        'charge_offset_x': 8,
        'charge_offset_y': -8,
        'subscript_offset_x': 4,
        'subscript_offset_y': 4,
    }
    
    # Style-specific settings
    style_specific = {
        'modern': {
            'atom_font_size': 16,
            'charge_font_size': 12,
            'subscript_font_size': 12,
            'bond_stroke_width': 2,
            'double_bond_spacing': 4,
            'triple_bond_spacing': 6,
            'min_bond_length': 20,
            'max_bond_length': 40,
        },
        'classic': {
            'atom_font_size': 14,
            'charge_font_size': 10,
            'subscript_font_size': 10,
            'bond_stroke_width': 1.5,
            'double_bond_spacing': 3,
            'triple_bond_spacing': 5,
            'min_bond_length': 25,
            'max_bond_length': 50,
        },
        'compact': {
            'atom_font_size': 12,
            'charge_font_size': 9,
            'subscript_font_size': 9,
            'bond_stroke_width': 1,
            'double_bond_spacing': 3,
            'triple_bond_spacing': 4,
            'min_bond_length': 15,
            'max_bond_length': 30,
        }
    }
    
    # If requested style doesn't exist, fall back to modern
    if style_name not in style_specific:
        style_name = 'modern'
    
    # Merge base settings with style-specific settings
    settings = {**base_settings, **style_specific[style_name]}
    
    return settings

def get_element_template(element_type: str, style_name: str = 'modern') -> str:
    """
    Get an SVG template for a specific molecular element.
    
    Args:
        element_type: Type of element ('atom', 'bond', 'charge', etc.)
        style_name: Name of the style ('modern', 'classic')
        
    Returns:
        SVG template string for the specified element
    """
    # Get style settings
    style = get_style_settings(style_name)
    
    # Templates for different element types
    templates = {
        'atom': f'''
            <text 
                x="{{}}" 
                y="{{}}" 
                font-family="{style['font_family']}" 
                font-size="{style['atom_font_size']}" 
                text-anchor="middle" 
                dominant-baseline="central" 
                fill="{{}}">{{}}</text>
        ''',
        
        'single_bond': f'''
            <line 
                x1="{{}}" 
                y1="{{}}" 
                x2="{{}}" 
                y2="{{}}" 
                stroke="{style['bond_color']}" 
                stroke-width="{style['bond_stroke_width']}" />
        ''',
        
        'double_bond': f'''
            <line 
                x1="{{}}" 
                y1="{{}}" 
                x2="{{}}" 
                y2="{{}}" 
                stroke="{style['bond_color']}" 
                stroke-width="{style['bond_stroke_width']}" />
            <line 
                x1="{{}}" 
                y1="{{}}" 
                x2="{{}}" 
                y2="{{}}" 
                stroke="{style['bond_color']}" 
                stroke-width="{style['bond_stroke_width']}" />
        ''',
        
        'triple_bond': f'''
            <line 
                x1="{{}}" 
                y1="{{}}" 
                x2="{{}}" 
                y2="{{}}" 
                stroke="{style['bond_color']}" 
                stroke-width="{style['bond_stroke_width']}" />
            <line 
                x1="{{}}" 
                y1="{{}}" 
                x2="{{}}" 
                y2="{{}}" 
                stroke="{style['bond_color']}" 
                stroke-width="{style['bond_stroke_width']}" />
            <line 
                x1="{{}}" 
                y1="{{}}" 
                x2="{{}}" 
                y2="{{}}" 
                stroke="{style['bond_color']}" 
                stroke-width="{style['bond_stroke_width']}" />
        ''',
        
        'charge': f'''
            <text 
                x="{{}}" 
                y="{{}}" 
                font-family="{style['font_family']}" 
                font-size="{style['charge_font_size']}" 
                text-anchor="middle" 
                fill="{{}}">{{}}</text>
        ''',
        
        'subscript': f'''
            <text 
                x="{{}}" 
                y="{{}}" 
                font-family="{style['font_family']}" 
                font-size="{style['subscript_font_size']}" 
                text-anchor="middle" 
                fill="{{}}">{{}}</text>
        ''',
        
        'benzene_ring': f'''
            <polygon 
                points="{{}}" 
                fill="none" 
                stroke="{style['bond_color']}" 
                stroke-width="{style['bond_stroke_width']}" />
            <circle 
                cx="{{}}" 
                cy="{{}}" 
                r="{{}}" 
                fill="none" 
                stroke="{style['bond_color']}" 
                stroke-width="1" 
                stroke-dasharray="2,2" />
        ''',
    }
    
    return templates.get(element_type, '')

def color_by_atom_type(element: str) -> str:
    """
    Get the color for a specific atom type.
    For now we're using black for all atoms per the specification,
    but this function allows for future customization.
    
    Args:
        element: Chemical element symbol
        
    Returns:
        Color as a string (hex or name)
    """
    # Per specification, we're using black for all atoms
    # This function is included for potential future customization
    return 'black'
