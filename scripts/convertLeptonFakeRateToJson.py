#!/usr/bin/env python3
"""
Convert lepton fake rate ROOT files to correctionlib JSON format

Features:
- Simplified arguments: just -e ERA -l LEPTON
- Integrates both main and QCD fake rate files into single JSON
- Applies lepton-specific edge restrictions (electron: 15-50 GeV, muon: 10-50 GeV)
- Sets proper axis names (absScEta/absEta + ptCorr)
- Configures appropriate flow behavior (error for eta, clamp for pt overflow)
"""

import os
import re
import argparse
import uproot
import numpy as np
import correctionlib
import correctionlib.schemav2 as cs

def _read_hist2d(file_path, key_name):
    with uproot.open(file_path) as f:
        obj = f[key_name]
        # uproot TH2 provides to_numpy() returning (values, xedges, yedges)
        vals, xedges, yedges = obj.to_numpy()
        return vals, xedges.tolist(), yedges.tolist()


def _get_axis_description(axis_name, lepton_type):
    """Get appropriate description for axis based on name and lepton type"""
    if axis_name == "absScEta":
        return "Absolute supercluster pseudorapidity"
    elif axis_name == "scEta":
        return "Supercluster pseudorapidity"
    elif axis_name == "absEta":
        return "Absolute pseudorapidity"
    elif axis_name == "eta":
        return "Pseudorapidity"
    else:
        return axis_name

def _restrict_pt(values, eta_edges, pt_edges, pt_min, pt_max):
    # Find pt bin indices within range (inclusive)
    start = 0
    while start < len(pt_edges)-1 and pt_edges[start] < pt_min:
        start += 1
    end = len(pt_edges)-1
    while end > 0 and pt_edges[end] > pt_max:
        end -= 1
    # Ensure boundaries
    new_edges = pt_edges[start:end+1]
    if new_edges[0] > pt_min:
        new_edges = [pt_min] + new_edges
    if new_edges[-1] < pt_max:
        new_edges = new_edges + [pt_max]
    # Slice values correspondingly
    # values shape: [n_eta_bins, n_pt_bins] assuming uproot order (x=eta, y=pt)
    n_eta_bins = len(eta_edges)-1
    n_pt_bins_orig = len(pt_edges)-1
    s = start
    e = min(end, n_pt_bins_orig-1)
    sliced = values[:, s:e]
    # If we inserted new boundary at start or end, pad by edge rows
    n_pt_bins_new = len(new_edges)-1
    if sliced.shape[1] < n_pt_bins_new:
        pad_left = 1 if new_edges[0] == pt_min and (start == 0 or pt_edges[start] > pt_min) else 0
        pad_right = 1 if new_edges[-1] == pt_max and (end == len(pt_edges)-1 or pt_edges[end] < pt_max) else 0
        if pad_left:
            sliced = np.concatenate([sliced[:, :1], sliced], axis=1)
        if pad_right and sliced.shape[1] < n_pt_bins_new:
            sliced = np.concatenate([sliced, sliced[:, -1:]], axis=1)
    return sliced, new_edges


def _make_correction(name, desc, era, lepton_type, values, eta_edges, pt_edges):
    # Use eta for Run3, abseta for Run2
    run3_eras = ['2022', '2022EE', '2023', '2023BPix']
    is_run3 = era in run3_eras

    #if is_run3:
    #    axis0 = "scEta" if lepton_type == "electron" else "eta"
    #else:
    #    axis0 = "absScEta" if lepton_type == "electron" else "absEta"
    axis0 = "absScEta" if lepton_type == "electron" else "absEta"
    axis1 = "ptCorr"
    # Flatten in eta-major, then pt bins
    n_eta = len(eta_edges)-1
    n_pt = len(pt_edges)-1
    content = []
    for i in range(n_eta):
        for j in range(n_pt):
            content.append(float(values[i, j]))
    data = cs.MultiBinning(
        nodetype="multibinning",
        inputs=[axis0, axis1],
        edges=[eta_edges, pt_edges],
        content=content,
        flow="clamp",
    )
    corr = cs.Correction(
        name=name,
        version=1,
        description=desc,
        inputs=[
            cs.Variable(name=axis0, type="real", description=_get_axis_description(axis0, lepton_type)),
            cs.Variable(name=axis1, type="real", description="Corrected transverse momentum [GeV]"),
        ],
        output=cs.Variable(name="out", type="real", description="out"),
        data=data,
    )
    return corr


def convert_fakerate_to_json(era, lepton_type):
    """Convert fake rate ROOT files to correctionlib JSON format with simplified arguments
    
    Args:
        era: Data-taking era (e.g., 2017, 2018, 2022)
        lepton_type: Lepton type (electron or muon)
    
    Returns:
        tuple: (output_file_path, number_of_corrections)
    """
    
    # Define file paths based on era and lepton type
    if lepton_type == "electron":
        lepton_dir = "EGM"
    else:
        lepton_dir = "MUO"
    
    base_path = f"data/Run3_v13_Run2_v9/{era}/{lepton_dir}/root"
    main_file = f"{base_path}/fakerate_TopHNT.root"
    mc_file = f"{base_path}/fakerate_MC_TopHNT.root"
    output_file = f"data/Run3_v13_Run2_v9/{era}/{lepton_dir}/fakerate_TopHNT.json"
    
    # Check if files exist
    if not os.path.exists(main_file):
        raise RuntimeError(f"Main fake rate file not found: {main_file}")
    if not os.path.exists(mc_file):
        raise RuntimeError(f"MC fake rate file not found: {mc_file}")
    
    corrections = []
    
    # Process main fake rate (Central)
    print(f"Processing main fake rate file: {main_file}")
    with uproot.open(main_file) as f:
        keys = [k.replace(';1','') for k in f.keys()]
    # Prefer key containing 'Central'
    key_central = next((k for k in keys if 'Central' in k or 'central' in k), keys[0])
    vals, eta_edges, pt_edges = _read_hist2d(main_file, key_central)
    # Assume uproot returns shape (n_eta, n_pt). If not, try transpose to match
    import numpy as np
    if vals.shape[0] != len(eta_edges)-1 or vals.shape[1] != len(pt_edges)-1:
        vals = vals.T
    pt_min, pt_max = (15.0, 100.0) if lepton_type == 'electron' else (10.0, 100.0)
    vals_r, pt_edges_r = _restrict_pt(vals, eta_edges, pt_edges, pt_min, pt_max)
    corr_central = _make_correction(
        name=f"fakerate_{lepton_type}_Central",
        desc=f"{lepton_type} fake rate - Central",
        era=era,
        lepton_type=lepton_type,
        values=vals_r,
        eta_edges=eta_edges,
        pt_edges=pt_edges_r,
    )
    corrections.append(corr_central)
    
    # Process MC fake rate file (all histograms)
    print(f"Processing MC fake rate file: {mc_file}")
    with uproot.open(mc_file) as f:
        mc_keys = [k.replace(';1','') for k in f.keys()]

    for key in mc_keys:
        # Extract variation name from key like "fake rate - (QCD_EMEnriched_bjet)"
        # Pattern: "fake rate - (NAME)" -> NAME
        match = re.search(r'\(([^)]+)\)', key)
        if not match:
            print(f"  Skipping unrecognized key: {key}")
            continue
        variation_name = match.group(1)

        vals_mc, eta_edges_mc, pt_edges_mc = _read_hist2d(mc_file, key)
        if vals_mc.shape[0] != len(eta_edges_mc)-1 or vals_mc.shape[1] != len(pt_edges_mc)-1:
            vals_mc = vals_mc.T
        vals_mc_r, pt_edges_mc_r = _restrict_pt(vals_mc, eta_edges_mc, pt_edges_mc, pt_min, pt_max)

        corr_mc = _make_correction(
            name=f"fakerate_{lepton_type}_{variation_name}",
            desc=f"{lepton_type} fake rate - {variation_name}",
            era=era,
            lepton_type=lepton_type,
            values=vals_mc_r,
            eta_edges=eta_edges_mc,
            pt_edges=pt_edges_mc_r,
        )
        corrections.append(corr_mc)
        print(f"  Added: fakerate_{lepton_type}_{variation_name}")
    
    if not corrections:
        raise RuntimeError("No histograms were successfully converted")
    
    # Create correction set
    cset = correctionlib.schemav2.CorrectionSet(
    schema_version=2,
        description=f"{lepton_type} fake rate corrections for {era} (includes MC systematic sources)",
        corrections=corrections
    )
    
    # Write to JSON file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as fout:
        fout.write(cset.model_dump_json(exclude_unset=True, indent=2))
    
    print(f"Successfully wrote {len(corrections)} corrections to {output_file}")
    return output_file, len(corrections)

def main():
    """Main function with simplified argument parsing"""
    parser = argparse.ArgumentParser(
        description="Convert lepton fake rate ROOT files to JSON format",
        epilog="""
Examples:
  # Convert muon fake rates for 2017
  python3 convertLeptonFakeRateToJson.py -e 2017 -l muon
  
  # Convert electron fake rates for 2022
  python3 convertLeptonFakeRateToJson.py -e 2022 -l electron
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-e", "--era", required=True, 
                       help="Data-taking era (e.g., 2017, 2018, 2022)")
    parser.add_argument("-l", "--lepton", choices=["electron", "muon"], required=True, 
                       help="Lepton type (electron or muon)")
    
    args = parser.parse_args()
    
    try:
        output_file, num_corrections = convert_fakerate_to_json(args.era, args.lepton)
        print(f"\nConversion completed successfully!")
        print(f"Output: {output_file}")
        print(f"Total corrections: {num_corrections}")
        
        # Show what was converted
        run3_eras = ['2022', '2022EE', '2023', '2023BPix']
        is_run3 = args.era in run3_eras
        #if is_run3:
        #    axis_name = "scEta" if args.lepton == "electron" else "eta"
        #else:
        #    axis_name = "absScEta" if args.lepton == "electron" else "absEta"
        axis_name = "absScEta" if args.lepton == "electron" else "absEta"

        print(f"\nFeatures applied:")
        print(f"  Era: {args.era} ({'Run3' if is_run3 else 'Run2'})")
        print(f"  Axis names: {axis_name} + ptCorr")
        print(f"  Flow: clamp for both axes (eta, pt)")
        print(f"  Integrated: Main + MC systematic corrections (flavor-tagged)")
        
        return 0
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())
