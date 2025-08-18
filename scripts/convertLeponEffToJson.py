#!/usr/bin/env python3
import os
import argparse
import correctionlib.schemav2 as cs
import json
import rich
from ROOT import TFile

parser = argparse.ArgumentParser(description="Convert a ROOT file to JSON format.")
parser.add_argument("--inputFile", type=str, help="Input ROOT file path")
parser.add_argument("--inputHist", type=str, help="Input histogram name (base name for nominal). Can parse multiple histograms e.g. eff_data:eff_mc")
parser.add_argument("--csetName", type=str, default="", help="Name of the correction set")
parser.add_argument("--xaxis", type=str, default="abseta", help="Name of the x-axis variable")
parser.add_argument("--yaxis", type=str, default="pt", help="Name of the y-axis variable")
parser.add_argument("--outputName", type=str, default="", help="Name of the output values")
parser.add_argument("--categoryName", type=str, default="", help="Name of the category")
parser.add_argument("--outputPath", type=str, help="Output JSON file path")
parser.add_argument("--name", type=str, default="", help="Name of the correction")
parser.add_argument("--version", type=int, default=1, help="Version of the correction")
args = parser.parse_args()

def create_binning(hist):
    """Create a binning structure from a ROOT histogram"""
    # Get x-axis (abseta) edges
    edges0 = []
    for i in range(hist.GetNbinsX() + 1):
        edges0.append(hist.GetXaxis().GetBinLowEdge(i + 1))
    
    # Get y-axis (pt) edges
    edges1 = []
    for i in range(hist.GetNbinsY() + 1):
        edges1.append(hist.GetYaxis().GetBinLowEdge(i + 1))
    
    # Create content structure for each abseta bin
    content = []
    for i in range(hist.GetNbinsX()):
        # Get values for this abseta bin
        values = []
        for j in range(hist.GetNbinsY()):
            values.append(hist.GetBinContent(i + 1, j + 1))
        
        # Create binning for pt
        pt_binning = cs.Binning(
            nodetype="binning",
            input=args.yaxis,
            edges=edges1,
            content=values,
            flow="clamp"
        )
        content.append(pt_binning)
    
    # Create main binning for abseta
    return cs.Binning(
        nodetype="binning",
        input=args.xaxis,
        edges=edges0,
        content=content,
        flow="clamp"
    )

# Open the ROOT file and get the histograms
with TFile.Open(args.inputFile) as f:
    # split inputHist by :
    inputHists = args.inputHist.split(":")
    outputNames = args.outputName.split(":")
    cset_items = []
    for inputHist, outputName in zip(inputHists, outputNames):
        # Get nominal histogram
        hist_nom = f.Get(inputHist)
        
        # Create up/down variations using bin errors instead of separate histograms
        hist_up = hist_nom.Clone(f"{inputHist}_up")
        hist_down = hist_nom.Clone(f"{inputHist}_down")
        
        # Fill up/down histograms by adding/subtracting bin errors
        for i in range(1, hist_nom.GetNbinsX() + 1):
            for j in range(1, hist_nom.GetNbinsY() + 1):
                nominal_value = hist_nom.GetBinContent(i, j)
                bin_error = hist_nom.GetBinError(i, j)
                hist_up.SetBinContent(i, j, nominal_value + bin_error)
                hist_down.SetBinContent(i, j, nominal_value - bin_error)
    
        # Create binning structures
        binning_nom = create_binning(hist_nom)
        binning_up = create_binning(hist_up)
        binning_down = create_binning(hist_down)

        cset = cs.Correction(
            name=inputHist,
            version=args.version,
            description=args.categoryName,
            inputs=[
                cs.Variable(name=args.xaxis, type="real", description=args.xaxis),
                cs.Variable(name=args.yaxis, type="real", description=args.yaxis),
                cs.Variable(name="syst", type="string", description="systematic variation (nom, up, down)")
            ],
            output=cs.Variable(name=outputName, type="real", description=outputName),
            data=cs.Category(
                nodetype="category",
                input="syst",
                content=[
                    cs.CategoryItem(key="nom", value=binning_nom),
                    cs.CategoryItem(key="up", value=binning_up),
                    cs.CategoryItem(key="down", value=binning_down)
                ]
            )
        )
        cset_items.append(cset)

# Create a correction set wrapper
corrections = cs.CorrectionSet(
    schema_version=2,
    description=args.csetName,
    corrections=cset_items
)

# Print the correction set
rich.print(corrections)

# Write to output file if specified
if args.outputPath:
    with open(args.outputPath, 'w') as f:
        json.dump(corrections.model_dump(), f, indent=2)
