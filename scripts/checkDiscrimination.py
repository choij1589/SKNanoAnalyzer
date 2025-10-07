#!/usr/bin/env python3
import sys
import ROOT

def checkDiscrimination(histfile, channel="SR3Mu"):
    """
    Extract discrimination accuracies from SignalKinematics output

    Histogram structure:
    - Bin 0: Incorrect selections
    - Bin 1: Correct selections

    Accuracy = Bin1 / (Bin0 + Bin1)
    """

    f = ROOT.TFile.Open(histfile)
    if not f or f.IsZombie():
        print(f"Error: Cannot open {histfile}")
        return

    variables = [
        "acoplanarity",
        "scalarPtSum",
        "ptAsymmetry",
        "gammaFactor_smaller",
        "gammaFactor_larger",
        "gammaAcop_smaller",
        "gammaAcop_larger",
        "deltaR_pair_mu3rd_larger",
        "deltaR_pair_mu3rd_smaller",
        "deltaPhi_pair_mu3rd_larger",
        "ptRatio_mu3rd_smaller",
        "deltaR_nearestBjet_smaller",
        "deltaR_nearestBjet_larger",
        "deltaPhi_pair_MET_smaller",
        "deltaPhi_pair_MET_larger",
        "MT_pair_MET_smaller",
        "MT_pair_MET_larger",
        "deltaR_leadingNonBjet_smaller",
        "deltaR_leadingNonBjet_larger",
        "deltaPhi_muSS_MET_larger",
        "deltaPhi_muSS_MET_smaller",
        "MT_muSS_MET_larger",
        "MT_muSS_MET_smaller",
        "MT_asymmetry_smaller",
        "MT_asymmetry_larger"
    ]

    print(f"\n{'='*60}")
    print(f"Discrimination Power Analysis for {channel}")
    print(f"File: {histfile}")
    print(f"{'='*60}\n")

    print(f"{'Variable':<20} {'Correct':<12} {'Total':<12} {'Accuracy':<12}")
    print(f"{'-'*60}")

    results = {}
    for var in variables:
        histpath = f"{channel}/Discrimination/{var}_correct"
        h = f.Get(histpath)

        if not h:
            print(f"{var:<20} {'NOT FOUND':<12}")
            continue

        # Get bin contents
        # Bin 0 = incorrect (x value 0-1)
        # Bin 1 = correct (x value 1-2)
        correct = h.GetBinContent(2)  # bin 1 (x=1)
        incorrect = h.GetBinContent(1)  # bin 0 (x=0)
        total = correct + incorrect

        if total > 0:
            accuracy = correct / total
            results[var] = accuracy
            print(f"{var:<20} {correct:<12.1f} {total:<12.1f} {accuracy:<12.3%}")
        else:
            print(f"{var:<20} {'NO DATA':<12}")

    if results:
        print(f"\n{'-'*60}")
        best_var = max(results, key=results.get)
        worst_var = min(results, key=results.get)
        print(f"\nBest discriminator:  {best_var:<30} ({results[best_var]:.1%})")
        print(f"Worst discriminator: {worst_var:<30} ({results[worst_var]:.1%})")
        print(f"\nInterpretation:")

        # Determine selection direction
        if "larger" in best_var:
            direction = "LARGER"
            var_name = best_var.replace("_larger", "")
        elif "smaller" in best_var:
            direction = "SMALLER"
            var_name = best_var.replace("_smaller", "")
        else:
            direction = "SMALLER"
            var_name = best_var

        print(f"  - Pick pair with {direction} {var_name}")
        print(f"  - This correctly identifies signal pair {results[best_var]:.1%} of the time")

    print(f"\n{'='*60}\n")

    f.Close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python checkDiscrimination.py <histfile.root> [channel]")
        print("Example: python checkDiscrimination.py hist_0.root SR3Mu")
        sys.exit(1)

    histfile = sys.argv[1]
    channel = sys.argv[2] if len(sys.argv) > 2 else "SR3Mu"

    checkDiscrimination(histfile, channel)
