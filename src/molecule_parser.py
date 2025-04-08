"""
Molecule Parser Module
Handles parsing of molecular input (formulas/SMILES)
"""

from typing import Optional, Dict, Any, Tuple, List
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def parse_molecule(molecule_input: str) -> Optional[Chem.Mol]:
    """
    Parse a molecule from various input formats.
    
    Args:
        molecule_input: SMILES string, molecular formula, or common name
        
    Returns:
        RDKit molecule object if parsing successful, None otherwise
    """
    # Try parsing as SMILES
    molecule = parse_smiles(molecule_input)
    if molecule:
        return molecule
    
    # If not SMILES, attempt to parse as a name using chemical name resolvers
    # (This would require additional resources - for now, let's assume SMILES only)
    
    return None

def parse_smiles(smiles: str) -> Optional[Chem.Mol]:
    """
    Parse a SMILES string into an RDKit molecule.
    
    Args:
        smiles: SMILES notation of a molecule
        
    Returns:
        RDKit molecule object if parsing successful, None otherwise
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if not molecule:
            return None
        
        # Generate 2D coordinates for visualization
        AllChem.Compute2DCoords(molecule)
        
        return molecule
    except:
        return None

def get_molecule_data(molecule: Chem.Mol) -> Dict[str, Any]:
    """
    Extract relevant data from a molecule for SVG generation.
    
    Args:
        molecule: RDKit molecule object
        
    Returns:
        Dictionary with atom coordinates, bond information, etc.
    """
    # Initialize data structure
    data = {
        'atoms': [],
        'bonds': [],
        'dimensions': {'width': 0, 'height': 0}
    }
    
    # Get conformer (2D coordinates)
    conf = molecule.GetConformer()
    
    # Process atoms
    for atom in molecule.GetAtoms():
        atom_idx = atom.GetIdx()
        pos = conf.GetAtomPosition(atom_idx)
        
        atom_data = {
            'idx': atom_idx,
            'symbol': atom.GetSymbol(),
            'x': pos.x,
            'y': pos.y,
            'charge': atom.GetFormalCharge(),
            'is_aromatic': atom.GetIsAromatic(),
            'num_h': atom.GetTotalNumHs()
        }
        
        data['atoms'].append(atom_data)
    
    # Process bonds
    for bond in molecule.GetBonds():
        bond_data = {
            'begin_atom_idx': bond.GetBeginAtomIdx(),
            'end_atom_idx': bond.GetEndAtomIdx(),
            'bond_type': bond.GetBondType(),
            'is_aromatic': bond.GetIsAromatic(),
            'is_conjugated': bond.GetIsConjugated(),
            'is_in_ring': bond.IsInRing()
        }
        
        data['bonds'].append(bond_data)
    
    # Calculate bounding box for the molecule
    x_coords = [atom['x'] for atom in data['atoms']]
    y_coords = [atom['y'] for atom in data['atoms']]
    
    if x_coords and y_coords:
        min_x, max_x = min(x_coords), max(x_coords)
        min_y, max_y = min(y_coords), max(y_coords)
        
        # Add some padding
        padding = 0.5
        width = max_x - min_x + padding * 2
        height = max_y - min_y + padding * 2
        
        data['dimensions']['width'] = width
        data['dimensions']['height'] = height
        data['dimensions']['min_x'] = min_x - padding
        data['dimensions']['min_y'] = min_y - padding
    
    return data

def parse_name_to_smiles(name: str) -> Optional[str]:
    """
    Convert a common chemical name to SMILES.
    This is a placeholder - would need to integrate with a chemical name database.
    
    Args:
        name: Common name of a chemical compound
        
    Returns:
        SMILES string if conversion successful, None otherwise
    """
    # This would ideally connect to PubChem, ChemSpider, or similar database
    # For now, we'll implement a small dictionary of common molecules
    common_molecules = {
        'methanol': 'CO',
        'ethanol': 'CCO',
        'water': 'O',
        'benzene': 'c1ccccc1',
        'acetic acid': 'CC(=O)O',
        'acetone': 'CC(=O)C',
        # Add more as needed
    }
    
    return common_molecules.get(name.lower())
