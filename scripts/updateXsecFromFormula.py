#!/usr/bin/env python3
import os
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--era', type=str, required=True, help='era (e.g., 2017)')
parser.add_argument('--dry-run', action='store_true', help='Show what would be updated without modifying files')
args = parser.parse_args()

def main():
    sampleInfoJson = os.path.join(os.environ['SKNANO_DATA'], args.era, 'Sample', 'CommonSampleInfo.json')

    # Load the sample information
    with open(sampleInfoJson, 'r') as f:
        sampleInfos = json.load(f)

    updated_count = 0

    # Process each sample
    for alias, sampleInfo in sampleInfos.items():
        # Only process MC samples with xsec_formula
        if sampleInfo.get("isMC") == 1 and "xsec_formula" in sampleInfo:
            old_xsec = sampleInfo.get("xsec", "N/A")
            formula = sampleInfo["xsec_formula"]

            try:
                # Evaluate the formula
                new_xsec = eval(formula)

                # Check if update is needed
                if old_xsec != new_xsec:
                    print(f"{alias}:")
                    print(f"  Formula: {formula}")
                    print(f"  Old xsec: {old_xsec}")
                    print(f"  New xsec: {new_xsec}")

                    if not args.dry_run:
                        sampleInfo["xsec"] = new_xsec

                    updated_count += 1
            except Exception as e:
                print(f"ERROR evaluating formula for {alias}: {formula}")
                print(f"  Error: {e}")

    # Save the updated file if not dry-run
    if not args.dry_run and updated_count > 0:
        with open(sampleInfoJson, 'w') as f:
            json.dump(sampleInfos, f, indent=4)
        print(f"\n✓ Updated {updated_count} samples in {sampleInfoJson}")
    elif args.dry_run:
        print(f"\n[DRY RUN] Would update {updated_count} samples")
    else:
        print(f"\nNo updates needed")

if __name__ == "__main__":
    main()
