from DataLoader import *
from hmm.Model import *
from hmm.ModelParams import *
import numpy as np

if __name__ == '__main__' :
    datafile = '8_DATA.PCA_normalized.filtered.sample_zscores.RD.txt'
    outputfile = 'output.txt'
    paramsfile = 'params.txt'

    dataloader = DataLoader(datafile)
    params = dataloader.getParams(paramsfile)
    dataloader.skipHeadline()
    sample = dataloader.getNextSample()
    while sample :
        targets = dataloader.getTargets()
        modelParams = ModelParams(params, targets)

        model = Model(modelParams, sample['observations'])
        pathlist = model.forwardBackward_Viterbi()
        dataloader.outputCNV(sample['sample_id'], targets, pathlist)
        sample = dataloader.getNextSample()