from DataLoader import *
from hmm.Model import *
from hmm.ModelParams import *
import numpy as np

if __name__ == '__main__' :
    datafile = '8_DATA.PCA_normalized.filtered.sample_zscores.RD.txt'
    #datafile = 'data'
    outputfile = 'output'
    paramsfile = 'params.txt'

    print 'Loading data file...'
    dataloader = DataLoader(datafile)
    params = dataloader.getParams(paramsfile)
    dataloader.skipHeadline()
    sample = dataloader.getNextSample()
    targets_list = dataloader.getTargetsList()
    output = file(outputfile, 'w')
    while sample :
        #target_index is used to split observations sequence
        target_index_begin = 0
        target_index_end = 0
        temp = 1
        for targets in targets_list:
            print 'Running HMM for sample[' + sample['sample_id'] + ']: ',
            print 'chr' + targets[0]._chr + ' [' + str(temp) + '\\' + str(len(targets_list)) + ']'
            temp += 1
            target_index_end = target_index_begin + len(targets)

            modelParams = ModelParams(params, targets)
            #the 'observations' of sample is splitted
            model = Model(modelParams, sample['observations'][target_index_begin:target_index_end])
            pathlist = model.forwardBackward_Viterbi()
            dataloader.outputCNV(output, sample['sample_id'], targets, pathlist, sample['observations'][target_index_begin:target_index_end])
            target_index_begin = target_index_end
        sample = dataloader.getNextSample()

    output.close()
