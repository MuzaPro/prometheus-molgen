"""
File Manager Module
Handles file operations, naming conventions, and directory management
"""

import os
import csv
import json
from typing import List, Tuple, Dict, Any

def save_svg(svg_content: str, name: str, output_dir: str = 'output') -> str:
    """
    Save SVG content to a file with proper naming convention.
    
    Args:
        svg_content: SVG content as a string
        name: Base name for the file (without extension)
        output_dir: Directory to save the file
        
    Returns:
        Path to the saved file
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Apply kebab-case convention
    filename = to_kebab_case(name) + '.svg'
    filepath = os.path.join(output_dir, filename)
    
    # Write content to file
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(svg_content)
    
    return filepath

def to_kebab_case(name: str) -> str:
    """
    Convert a string to kebab-case.
    
    Args:
        name: Input string
        
    Returns:
        Kebab-case string
    """
    # Replace spaces and special characters with hyphens
    name = name.lower().strip()
    # Replace multiple spaces and special chars with a single hyphen
    import re
    name = re.sub(r'[^a-z0-9]+', '-', name)
    # Remove leading/trailing hyphens
    name = name.strip('-')
    return name

def read_molecule_list(file_path: str) -> List[Tuple[str, str]]:
    """
    Read a list of molecules from a file.
    Supported formats: CSV, JSON, TXT
    
    Args:
        file_path: Path to the input file
        
    Returns:
        List of tuples (name, smiles)
    """
    ext = os.path.splitext(file_path)[1].lower()
    
    if ext == '.csv':
        return read_csv_molecules(file_path)
    elif ext == '.json':
        return read_json_molecules(file_path)
    elif ext == '.txt':
        return read_txt_molecules(file_path)
    else:
        raise ValueError(f"Unsupported file format: {ext}")

def read_csv_molecules(file_path: str) -> List[Tuple[str, str]]:
    """
    Read molecules from a CSV file.
    Expected format: name,smiles
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        List of tuples (name, smiles)
    """
    molecules = []
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        # Check for header
        first_row = next(reader)
        
        # If first row doesn't look like a SMILES string, treat it as header
        if len(first_row) >= 2 and not any(c in first_row[1] for c in '()[]'):
            # Skip, it was a header
            pass
        else:
            # Not a header, process the first row
            if len(first_row) >= 2:
                molecules.append((first_row[0], first_row[1]))
        
        # Process remaining rows
        for row in reader:
            if len(row) >= 2:
                molecules.append((row[0], row[1]))
    
    return molecules

def read_json_molecules(file_path: str) -> List[Tuple[str, str]]:
    """
    Read molecules from a JSON file.
    Expected format: [{"name": "...", "smiles": "..."}, ...]
    
    Args:
        file_path: Path to the JSON file
        
    Returns:
        List of tuples (name, smiles)
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    molecules = []
    
    # Handle different possible JSON structures
    if isinstance(data, list):
        for item in data:
            if isinstance(item, dict) and 'name' in item and 'smiles' in item:
                molecules.append((item['name'], item['smiles']))
    elif isinstance(data, dict):
        if 'molecules' in data and isinstance(data['molecules'], list):
            for item in data['molecules']:
                if isinstance(item, dict) and 'name' in item and 'smiles' in item:
                    molecules.append((item['name'], item['smiles']))
        else:
            # Treat as a dictionary of name: smiles pairs
            for name, smiles in data.items():
                if isinstance(smiles, str):
                    molecules.append((name, smiles))
    
    return molecules

def read_txt_molecules(file_path: str) -> List[Tuple[str, str]]:
    """
    Read molecules from a plain text file.
    Expected format: Each line contains "name smiles" or "name,smiles"
    
    Args:
        file_path: Path to the text file
        
    Returns:
        List of tuples (name, smiles)
    """
    molecules = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Check if comma-separated
            if ',' in line:
                parts = line.split(',', 1)
                if len(parts) >= 2:
                    name, smiles = parts[0].strip(), parts[1].strip()
                    molecules.append((name, smiles))
            else:
                # Assume space-separated
                parts = line.split(maxsplit=1)
                if len(parts) >= 2:
                    name, smiles = parts[0].strip(), parts[1].strip()
                    molecules.append((name, smiles))
    
    return molecules

def export_to_question_json(molecules: List[Dict[str, Any]], output_file: str) -> None:
    """
    Export molecule data to a format compatible with the quiz game's question.json.
    
    Args:
        molecules: List of molecule data dictionaries
        output_file: Path to output JSON file
    """
    # Convert molecule data to the format expected by the quiz game
    questions_data = {
        "questions": []
    }
    
    for i, mol in enumerate(molecules, 1):
        question = {
            "id": f"q{i}",
            "questionText": f"Identify the following molecule:",
            "product": {
                "name": mol["name"],
                "formula": mol.get("formula", ""),
                "imagePath": f"assets/molecules/{to_kebab_case(mol['name'])}.svg",
                "textRepresentation": mol.get("textRepresentation", mol.get("smiles", ""))
            },
            # Other question properties would be added here
        }
        questions_data["questions"].append(question)
    
    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(questions_data, f, indent=2)

def validate_svg_files(directory: str) -> Dict[str, List[str]]:
    """
    Validate SVG files against specifications.
    
    Args:
        directory: Directory containing SVG files
        
    Returns:
        Dictionary of validation results
    """
    results = {
        "valid": [],
        "invalid": [],
        "errors": {}
    }
    
    # List SVG files
    svg_files = [f for f in os.listdir(directory) if f.endswith('.svg')]
    
    for svg_file in svg_files:
        filepath = os.path.join(directory, svg_file)
        
        try:
            # Read file content
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Basic validation
            if not content.strip().startswith('<svg'):
                results["invalid"].append(svg_file)
                results["errors"][svg_file] = "Not a valid SVG file"
                continue
            
            # Check for required attributes
            if 'viewBox' not in content:
                results["invalid"].append(svg_file)
                results["errors"][svg_file] = "Missing viewBox attribute"
                continue
            
            # Add more validation checks as needed
            
            # If passes all checks
            results["valid"].append(svg_file)
            
        except Exception as e:
            results["invalid"].append(svg_file)
            results["errors"][svg_file] = str(e)
    
    return results
