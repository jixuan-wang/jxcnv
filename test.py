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
    targets_list = dataloader.getTargetsList()
    while sample :
        output = file(sample['sample_id'], 'w')
        #target_index is used to split observations sequence
        target_index_begin = 0
        target_index_end = 0
        for targets in targets_list:
            target_index_end = target_index_begin + len(targets)

            modelParams = ModelParams(params, targets)
            #the 'observations' of sample is splitted
            model = Model(modelParams, sample['observations'][target_index_begin:target_index_end])
            pathlist = model.forwardBackward_Viterbi()
            dataloader.outputCNV(output, targets, pathlist)
            sample = dataloader.getNextSample()

            target_index_begin = target_index_end

        output.close()
