#!/bin/bash
RUN=$1
CHANNEL=$2
MHc=$3

## Update DATASTREAM
DATASTREAM=""
MC_LIST=""
if [[ $RUN == "Run2" ]]; then
    if [[ $CHANNEL == "Run1E2Mu" || $CHANNEL == "Run2E1Mu" ]]; then
        DATASTREAM="Skim_TriLep_MuonEG"
    elif [[ $CHANNEL == "Run3Mu" ]]; then
        DATASTREAM="Skim_TriLep_DoubleMuon"
    else
        echo "Unknown channel: $CHANNEL"
        exit 1
    fi
    MC_LIST="SampleLists/Run2NanoV9"
elif [[ $RUN == "Run3" ]]; then
    if [[ $CHANNEL == "Run1E2Mu" || $CHANNEL == "Run2E1Mu" ]]; then
        DATASTREAM="Skim_TriLep_MuonEG"
    elif [[ $CHANNEL == "Run3Mu" ]]; then
        DATASTREAM="Skim_TriLep_DoubleMuon,Skim_TriLep_Muon,Skim_TriLep_Muon0,Skim_TriLep_Muon1"
    else
        echo "Unknown channel: $CHANNEL"
        exit 1
    fi
    MC_LIST="SampleLists/Run3NanoV13"
else
    echo "Unknown run: $RUN"
    exit 1
fi

PN_MODEL_DIR="${SKNANO_DATA}/All/Combined/Classifiers/ParticleNetMD"

get_pn_signal_samples() {
    local mhc="$1"
    local points=()

    if [[ -n "${mhc}" ]]; then
        shopt -s nullglob
        for model_dir in "${PN_MODEL_DIR}/${mhc}"_MA*; do
            [[ -f "${model_dir}/best_model/model.pt" ]] || continue
            [[ -f "${model_dir}/best_model/model_info.json" ]] || continue
            points+=("$(basename "${model_dir}")")
        done
        shopt -u nullglob
    else
        points=("MHc100_MA95" "MHc130_MA90" "MHc160_MA85")
    fi

    if [[ ${#points[@]} -eq 0 ]]; then
        echo "No ParticleNetMD models found for ${mhc} under ${PN_MODEL_DIR}" >&2
        exit 1
    fi

    local samples=()
    for point in "${points[@]}"; do
        samples+=("TTToHcToWAToMuMu-${point}")
    done

    local IFS=,
    echo "${samples[*]}"
}

PN_SIGNAL_SAMPLES=$(get_pn_signal_samples "${MHc}")
PN_FLAGS="${CHANNEL}"
if [[ -n "${MHc}" ]]; then
    PN_FLAGS="${CHANNEL},${MHc}"
fi

SKNano.py -a PromptAnalyzer -i ${DATASTREAM} -n 10 -r ${RUN} --userflags ${PN_FLAGS},NoHistMode --python --memory 4000
SKNano.py -a MatrixAnalyzer -i ${DATASTREAM} -n 10 -r ${RUN} --userflags ${PN_FLAGS},NoHistMode --python --memory 4000
SKNano.py -a PromptAnalyzer -i ${MC_LIST}/TriLepton.txt -n 40 -r ${RUN} --userflags ${PN_FLAGS},RunSyst,NoHistMode --python --memory 8000
SKNano.py -a PromptAnalyzer -i ${PN_SIGNAL_SAMPLES} -n 20 -r ${RUN} --userflags ${PN_FLAGS},RunSyst,RunTheoryUnc,NoHistMode --python --memory 8000
