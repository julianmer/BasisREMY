"""
Utility functions for finding and filtering MRS data files

Combines file format detection with Excel-based sequence mapping.
Provides dynamic file discovery based on:
- File format (extension) - detected dynamically
- Sequence type (from Excel metadata)
"""

import os
import sys
import pandas as pd
from collections import defaultdict

# Add project root to path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


# Configuration for MRS file formats and sequences
MRS_EXTENSIONS = {'.spar', '.7', '.dat', '.rda', '.ima', 'method', 'acqp', '.nii.gz', '.nii'}
SKIP_EXTENSIONS = {'.xlsx', '.md', '.txt', '.png', '.jpg', '.pdf', '.pyc', '.DS_Store'}
SKIP_FILES = {'T1.nii.gz'}

# Sequence detection patterns (centralized - easy to extend)
SEQUENCE_PATTERNS = {
    'PRESS': lambda p: 'PRESS' in p and 'MEGA' not in p,
    'STEAM': lambda p: 'STEAM' in p,
    'sLASER': lambda p: 'SLASER' in p or 'OSLASER' in p,
    'LASER': lambda p: 'LASER' in p and 'SLASER' not in p,
    # TODO: Add more sequences here as needed
    # 'MEGA-PRESS': lambda p: 'MEGAPRESS' in p or 'MEGA-PRESS' in p,
    # 'HERMES': lambda p: 'HERMES' in p,
    # 'HERCULES': lambda p: 'HERCULES' in p,
}


def _detect_mrs_format(filename):
    """
    Detect MRS file format from filename

    Returns:
        format_extension (str) or None
    """
    file_lower = filename.lower()

    # Check known extensions
    if file_lower.endswith('.spar'):
        return '.spar'
    elif filename.endswith('.7'):
        return '.7'
    elif file_lower.endswith('.dat'):
        return '.dat'
    elif file_lower.endswith('.rda'):
        return '.rda'
    elif file_lower.endswith('.ima'):
        return '.ima'
    elif file_lower.endswith('method') or file_lower == 'method':
        return 'method'
    elif file_lower.endswith('acqp') or file_lower == 'acqp':
        return 'acqp'
    elif file_lower.endswith('.nii.gz'):
        return '.nii.gz'
    elif file_lower.endswith('.nii'):
        return '.nii'

    return None


def _detect_sequence_from_filename(filename):
    """
    Detect sequence from filename (for BigGABA and similar datasets)

    Returns:
        sequence_name or None
    """
    filename_upper = filename.upper()

    for sequence_name, pattern_func in SEQUENCE_PATTERNS.items():
        if pattern_func(filename_upper):
            return sequence_name

    return None


def _build_folder_to_sequence_map(example_data_dir, sequence_map):
    """
    Build mapping from actual folder names to sequences

    Handles the mapping between Excel dataset names (Dataset_B1, Dataset_G1)
    and actual folder names (Dataset_00_Bruker_14T_STEAM_08)

    Returns:
        dict: {actual_folder_name: sequence}
    """
    folder_to_sequence = {}
    vendor_counters = {'B': 0, 'G': 0, 'P': 0, 'S': 0}

    remy_dir = os.path.join(example_data_dir, 'REMY_tests')
    if not os.path.exists(remy_dir):
        return folder_to_sequence

    for folder in sorted(os.listdir(remy_dir)):
        if folder.startswith('Dataset_'):
            # Extract vendor from folder name
            # Dataset_00_Bruker -> B, Dataset_02_GE -> G, etc.
            parts = folder.split('_')
            if len(parts) >= 3:
                vendor_name = parts[2]
                vendor_code = None

                if vendor_name.startswith('Bruker'):
                    vendor_code = 'B'
                elif vendor_name.startswith('GE'):
                    vendor_code = 'G'
                elif vendor_name.startswith('Philips'):
                    vendor_code = 'P'
                elif vendor_name.startswith('Siemens'):
                    vendor_code = 'S'

                if vendor_code:
                    vendor_counters[vendor_code] += 1
                    mapped_name = f'Dataset_{vendor_code}{vendor_counters[vendor_code]}'

                    # Look up sequence for this mapped name
                    if mapped_name in sequence_map:
                        folder_to_sequence[folder] = sequence_map[mapped_name]

    return folder_to_sequence


def _get_sequence_for_file(filepath, example_data_dir, folder_to_sequence):
    """
    Determine sequence for a given file

    Returns:
        sequence_name or 'Unknown'
    """
    rel_path = os.path.relpath(filepath, example_data_dir)

    # Check if in REMY_tests - use folder mapping
    if 'REMY_tests' in rel_path:
        parts = rel_path.split(os.sep)
        for part in parts:
            if part.startswith('Dataset_'):
                sequence = folder_to_sequence.get(part)
                if sequence:
                    return sequence
                break

    # Check BigGABA - use filename
    if 'BigGABA' in rel_path:
        filename = os.path.basename(filepath)
        sequence = _detect_sequence_from_filename(filename)
        if sequence:
            return sequence

    return 'Unknown'


def _detect_sequence_from_protocol(protocol):
    """
    Detect sequence type from protocol description

    Args:
        protocol: Protocol description string

    Returns:
        Sequence name or None
    """
    protocol_upper = protocol.upper()

    for sequence_name, pattern_func in SEQUENCE_PATTERNS.items():
        if pattern_func(protocol_upper):
            return sequence_name

    return None


def load_sequence_mapping(excel_path='example_data/nbm70039-sup-0001-supplementary_material.xlsx'):
    """
    Load sequence mapping from Excel file using Protocol Description column

    Returns:
        dict: {dataset_folder_name: sequence_type}
    """
    if not os.path.exists(excel_path):
        return {}

    try:
        df = pd.read_excel(excel_path)

        # First row is the header row
        df.columns = df.iloc[0]
        df = df[1:]  # Remove header row

        mapping = {}
        current_vendor = None
        dataset_counter = {}  # Track dataset numbers per vendor

        for idx, row in df.iterrows():
            dataset = str(row.iloc[0]).strip()

            # Check if this is a vendor header row
            if dataset in ['Bruker', 'GE', 'Philips', 'Siemens']:
                current_vendor = dataset
                if current_vendor not in dataset_counter:
                    dataset_counter[current_vendor] = 0
                continue

            # Skip non-dataset rows
            if not dataset.isdigit():
                continue

            # Get Protocol Description (column index 6)
            protocol = str(row.iloc[6]).strip() if pd.notna(row.iloc[6]) else ''

            if not protocol or protocol == 'nan':
                continue

            # Create dataset folder name matching filesystem
            # Format: Dataset_<vendor_initial><number>
            # e.g., Dataset_B1, Dataset_G1, Dataset_P1, Dataset_S1
            if current_vendor:
                dataset_counter[current_vendor] += 1
                vendor_code = current_vendor[0]  # B, G, P, or S
                dataset_num = dataset_counter[current_vendor]
                dataset_folder = f'Dataset_{vendor_code}{dataset_num}'
            else:
                dataset_folder = f'Dataset_{dataset}'

            # Detect sequence using centralized logic
            sequence = _detect_sequence_from_protocol(protocol)

            if sequence:
                mapping[dataset_folder] = sequence

        return mapping

    except Exception as e:
        print(f"Warning: Could not load Excel sequence mapping: {e}")
        return {}


def find_all_mrs_files(example_data_dir='example_data', group_by='format'):
    """
    Dynamically find ALL MRS data files in example_data

    Args:
        example_data_dir: Path to example_data directory
        group_by: How to group files - 'format' or 'sequence' or 'both'

    Returns:
        If group_by='format': {extension: [filepaths]}
        If group_by='sequence': {sequence_name: [filepaths]}
        If group_by='both': {(extension, sequence): [filepaths]}
    """
    # Collect all MRS files
    all_files = []
    for root, dirs, filenames in os.walk(example_data_dir):
        # Skip hidden directories and __pycache__
        dirs[:] = [d for d in dirs if not d.startswith('.') and d != '__pycache__']

        for file in filenames:
            # Skip non-MRS files
            if file.startswith('.'):
                continue
            if file in SKIP_FILES:
                continue
            if any(file.endswith(ext) for ext in SKIP_EXTENSIONS):
                continue

            # Detect format
            detected_ext = _detect_mrs_format(file)
            if detected_ext:
                filepath = os.path.join(root, file)
                all_files.append((filepath, detected_ext))

    # Group by format only
    if group_by == 'format':
        files_by_format = defaultdict(list)
        for filepath, ext in all_files:
            files_by_format[ext].append(filepath)

        # Sort each list
        for ext in files_by_format:
            files_by_format[ext] = sorted(files_by_format[ext])

        return dict(files_by_format)

    # Group by sequence or both - need sequence mapping
    sequence_map = load_sequence_mapping()
    folder_to_sequence = _build_folder_to_sequence_map(example_data_dir, sequence_map)

    # Group by sequence only
    if group_by == 'sequence':
        files_by_sequence = defaultdict(list)

        for filepath, ext in all_files:
            sequence = _get_sequence_for_file(filepath, example_data_dir, folder_to_sequence)
            files_by_sequence[sequence].append(filepath)

        # Sort each list
        for seq in files_by_sequence:
            files_by_sequence[seq] = sorted(files_by_sequence[seq])

        return dict(files_by_sequence)

    # Group by both format AND sequence
    elif group_by == 'both':
        files_by_both = defaultdict(list)

        for filepath, ext in all_files:
            sequence = _get_sequence_for_file(filepath, example_data_dir, folder_to_sequence)
            key = (ext, sequence)
            files_by_both[key].append(filepath)

        # Sort each list
        for key in files_by_both:
            files_by_both[key] = sorted(files_by_both[key])

        return dict(files_by_both)

    else:
        raise ValueError(f"Invalid group_by: {group_by}. Must be 'format', 'sequence', or 'both'")


def get_files_for_format(format_ext, example_data_dir='example_data'):
    """
    Get all files for a specific format

    Args:
        format_ext: File extension like '.spar', '.7', '.dat', 'method', etc.
        example_data_dir: Path to example_data directory

    Returns:
        List of filepaths with that extension
    """
    all_files = find_all_mrs_files(example_data_dir, group_by='format')
    return all_files.get(format_ext, [])


def get_files_for_sequence(sequence_name, example_data_dir='example_data'):
    """
    Get all files for a specific sequence type

    Args:
        sequence_name: Sequence name like 'PRESS', 'STEAM', 'sLASER', etc.
        example_data_dir: Path to example_data directory

    Returns:
        List of filepaths with that sequence
    """
    all_files = find_all_mrs_files(example_data_dir, group_by='sequence')
    return all_files.get(sequence_name, [])


def get_files_for_format_and_sequence(format_ext, sequence_name, example_data_dir='example_data'):
    """
    Get all files matching BOTH a format AND a sequence

    Args:
        format_ext: File extension like '.spar', '.7', etc.
        sequence_name: Sequence name like 'PRESS', 'STEAM', etc.
        example_data_dir: Path to example_data directory

    Returns:
        List of filepaths matching both criteria
    """
    all_files = find_all_mrs_files(example_data_dir, group_by='both')
    return all_files.get((format_ext, sequence_name), [])


if __name__ == '__main__':
    """Test the utility functions"""
    print("="*80)
    print("Testing MRS File Discovery Utilities")
    print("="*80)

    # Test format grouping
    print("\n1. Files grouped by FORMAT:")
    by_format = find_all_mrs_files(group_by='format')
    for ext, files in sorted(by_format.items()):
        print(f"  {ext:12s}: {len(files):4d} files")

    # Test sequence grouping
    print("\n2. Files grouped by SEQUENCE:")
    by_sequence = find_all_mrs_files(group_by='sequence')
    for seq, files in sorted(by_sequence.items()):
        print(f"  {seq:12s}: {len(files):4d} files")

    # Test combined grouping
    print("\n3. Files grouped by FORMAT + SEQUENCE:")
    by_both = find_all_mrs_files(group_by='both')
    for (ext, seq), files in sorted(by_both.items()):
        print(f"  {ext:12s} + {seq:12s}: {len(files):4d} files")

    print("\n" + "="*80)















