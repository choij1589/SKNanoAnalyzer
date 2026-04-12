#!/bin/bash
# Fix stuck skim postproc for TTLL systematic samples
# Moves Temp_Skim_TriLep_* -> Skim_TriLep_* and generates skim info JSONs

set -e

# Requires setup.sh to have been sourced beforehand
if [ -z "$SKNANO_PYTHON" ] || [ -z "$SKNANO_DATA" ]; then
    echo "Error: SKNANO_PYTHON or SKNANO_DATA not set. Run 'source setup.sh' first."
    exit 1
fi

SKIM_BASE="/gv0/DATA/SKNano/Run2NanoAODv9p1"
ERAS=("2016preVFP" "2016postVFP" "2017" "2018")

processed=0
skipped=0

for era in "${ERAS[@]}"; do
    skim_dir="$SKIM_BASE/$era/MC/Skim/choij"
    if [ ! -d "$skim_dir" ]; then
        echo "[$era] Skim directory not found, skipping"
        continue
    fi

    for temp_dir in "$skim_dir"/Temp_Skim_TriLep_TTLL_*; do
        [ -d "$temp_dir" ] || continue

        dirname=$(basename "$temp_dir")
        final_name="${dirname/Temp_/}"
        final_dir="$skim_dir/$final_name"

        # Safety: skip if final directory already exists
        if [ -d "$final_dir" ]; then
            echo "[$era] SKIP: $final_name already exists"
            skipped=$((skipped + 1))
            continue
        fi

        # Move temp -> final
        echo "[$era] mv $dirname -> $final_name"
        mv "$temp_dir" "$final_dir"

        # Extract original PD name: Skim_TriLep_TTLL_foo -> TTLL_foo
        origPD="${final_name#Skim_TriLep_}"

        # Generate skim info JSON
        echo "[$era] sampleManager.py --makeSkimTreeInfo for $origPD"
        cd "$SKNANO_PYTHON"
        python3 sampleManager.py --era "$era" --makeSkimTreeInfo \
            --skimTreeFolder "$skim_dir" \
            --skimTreeSuffix TriLep \
            --skimTreeOrigPD "$origPD"

        processed=$((processed + 1))
    done
done

echo ""
echo "Done: $processed processed, $skipped skipped"
