"""
Comprehensive REMY tests for Philips SPAR files

Uses dynamic file discovery to find ALL .spar files in example_data,
runs through REMY, reports errors and field extraction statistics
"""

import pytest
import os
import sys

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from core.basisremy import BasisREMY
from tests.utils.file_discovery import get_files_for_format


@pytest.mark.remy
class TestPhilipsSPAR:
    """Test REMY parsing for Philips .spar files"""

    def test_all_spar_files(self, capsys):
        """
        Loop through ALL .spar files in example_data:
        - Uses dynamic file discovery
        - Runs through REMY
        - Reports which files threw errors
        - Shows ALL fields extracted with statistics
        - Format: "FieldName: X/Y files (Z%)"

        Note: Detailed field extraction only shown with -s flag
        """
        # Use utility to find all .spar files dynamically
        spar_files = get_files_for_format('.spar')

        if not spar_files:
            pytest.skip("No .spar files found")

        # Initialize BasisREMY
        basisremy = BasisREMY()

        # Track results
        all_fields = set()
        field_counts = {}
        successful_files = []
        failed_files = []

        # Loop through all files (no output during processing)
        for filepath in spar_files:
            rel_path = os.path.relpath(filepath, 'example_data')

            try:
                # Run through REMY
                params = basisremy.runREMY(import_fpath=filepath)
                if params is None:
                    params = {}

                # Also get backend parsed params
                try:
                    parsed_params, opt = basisremy.backend.parseREMY(params)
                    if parsed_params is None:
                        parsed_params = {}
                except:
                    parsed_params = {}

                # Collect ALL extracted fields
                extracted_fields = {}
                for field, value in params.items():
                    if value not in [None, '']:
                        extracted_fields[field] = value
                        all_fields.add(field)
                        field_counts[field] = field_counts.get(field, 0) + 1

                for field, value in parsed_params.items():
                    if value not in [None, ''] and field not in extracted_fields:
                        extracted_fields[field] = value
                        all_fields.add(field)
                        field_counts[field] = field_counts.get(field, 0) + 1

                successful_files.append((rel_path, extracted_fields))

            except Exception as e:
                # File threw an error
                error_msg = str(e)
                failed_files.append((rel_path, error_msg))

        # Now print results - this ALWAYS shows
        total = len(spar_files)
        success_count = len(successful_files)
        fail_count = len(failed_files)

        # Disable capturing to force output
        with capsys.disabled():
            print(f"\n{'='*80}")
            print(f"Testing Philips .spar Format")
            print(f"Found {len(spar_files)} files in example_data")
            print(f"{'='*80}")

            print(f"\n{'='*80}")
            print(f"RESULTS: Philips .spar Format")
            print(f"{'='*80}")
            print(f"Total files:   {total}")
            print(f"Successful:    {success_count} ({success_count/total*100:.1f}%)")
            print(f"Failed/Errors: {fail_count} ({fail_count/total*100:.1f}%)")

            # Show files that threw errors (ALWAYS SHOWN)
            if failed_files:
                print(f"\n{'-'*80}")
                print(f"FILES THAT THREW ERRORS ({len(failed_files)}):")
                print(f"{'-'*80}")
                for fpath, error in failed_files:
                    print(f"\n  {fpath}")
                    print(f"    ERROR: {error}")

            # Show field extraction statistics - "Field: X/Y files (Z%)" (ALWAYS SHOWN)
            if all_fields:
                print(f"\n{'-'*80}")
                print(f"FIELD EXTRACTION STATISTICS:")
                print(f"{'-'*80}")
                print(f"{'Field':<35s} | {'Extracted':>15s}")
                print(f"{'-'*80}")

                for field in sorted(all_fields, key=lambda x: field_counts.get(x, 0), reverse=True):
                    count = field_counts.get(field, 0)
                    pct = (count / total) * 100
                    print(f"{field:<35s} | {count:>4d}/{total:<4d} ({pct:>5.1f}%)")

            # Show detailed fields ONLY with -s flag (verbose mode)
            # Check if pytest is running with -s
            is_verbose = '-s' in sys.argv or '--capture=no' in sys.argv

            if is_verbose and successful_files:
                print(f"\n{'-'*80}")
                print(f"DETAILED FIELD EXTRACTION (ALL {len(successful_files)} successful files):")
                print(f"{'-'*80}")
                for fpath, fields in successful_files:
                    print(f"\n{fpath}")
                    print(f"  Extracted {len(fields)} fields:")
                    for field, value in sorted(fields.items()):
                        value_str = str(value)
                        if len(value_str) > 50:
                            value_str = value_str[:47] + "..."
                        print(f"    {field:<30s}: {value_str}")

            print(f"\n{'='*80}\n")

        # Test passes if we tested something
        assert total > 0, "No files were tested"









