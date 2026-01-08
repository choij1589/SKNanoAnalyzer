#!/bin/bash
RUN=$1

SKNano.py -a EvtTreeProducer -i SampleLists/Run2NanoV9/TrainDataset.txt -r $RUN -n 10 --userflags Run1E2Mu --memory 4000
SKNano.py -a EvtTreeProducer -i SampleLists/Run2NanoV9/TrainDataset.txt -r $RUN -n 10 --userflags Run3Mu --memory 4000
