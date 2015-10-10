import argparse
from DataLoader import *
from hmm.Model import *
from hmm.ModelParams import *
import numpy as np

def svd(args):
	filename = args.datafile 
	
	# count the columns number of the data file
	f = open(filename)
	line = f.readline()
	colsnum = len(line.split('\t'))

	# skip 1st row and 1st column
	data = np.loadtxt(filename, dtype=np.float, delimiter='\t', skiprows=1, usecols=range(1, colsnum)) 
	# test for conifer
	# data = np.loadtxt(filename, dtype=np.float, delimiter='\t')

	### for test
	data = np.transpose(data) 

	# svd transform
	# comp_removed = args.svd
	U, S, Vt = np.linalg.svd(data, full_matrices=False)
	
	# new_S = np.diag(np.hstack([np.zeros([comp_removed]), S[comp_removed:]]))
	index = S < 0.7 * np.mean(S)
	new_S = np.diag(S * index)
	
	# reconstruct data matrix
	data = np.dot(U, np.dot(new_S, Vt))
	   
	# save to files
	file_u = open(args.output + '.U', 'w')
	file_s = open(args.output + '.S', 'w')
	file_vt = open(args.output + '.Vt', 'w')
	file_svd = open(args.output + '.SVD', 'w')

	np.savetxt(file_u, U, delimiter='\t')
	np.savetxt(file_s, S, delimiter='\t')
	np.savetxt(file_vt, Vt, delimiter='\t')
	np.savetxt(file_svd, np.transpose(data), delimiter='\t')


	file_u.close()
	file_s.close()
	file_vt.close()
	file_svd.close()	
	
def discover(args) :
    datafile = args.datafile
    outputfile = args.output
    paramsfile = args.params
    sample_req = args.sample
    sample_flag = False #used to check whether sample_req exists

    print 'Loading data file...'
    dataloader = DataLoader(datafile)
    params = dataloader.getParams(paramsfile)
    dataloader.skipHeadline()
    sample = dataloader.getNextSample()
    targets_list = dataloader.getTargetsList()
    output = file(outputfile, 'w')
    output.write('SAMPLE_ID\tCNV\tFULL_INTERVAL\tINDEX\tINTERVAL\tREAD_DEPTH\n')
    while sample :
        if sample_req == '' or (sample_req != '' and sample['sample_id'] == sample_req):
            sample_flag = True
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
    dataloader.closeFile()

    if not sample_flag:
        print 'Could not find the sample_id specified.'

parser = argparse.ArgumentParser(prog='jxcnv', description='Designed by jx.')
subparsers = parser.add_subparsers()

#SVD
svd_parser = subparsers.add_parser('svd', help="SVD")
svd_parser.add_argument('--datafile', required=True, help='')
svd_parser.add_argument('--output', required=True, help='')
# svd_parser.add_argument('--svd', type=int, required=True, help='Number of components to remove')
svd_parser.set_defaults(func=svd)

#CNV discover
cnv_parser = subparsers.add_parser('discover', help="Run HMM to discover CNVs")
cnv_parser.add_argument('--params', required=True, help='Parameters used by HMM')
cnv_parser.add_argument('--datafile', required=True, help='Read depth file.')
cnv_parser.add_argument('--output', required=True, help='Output file.')
cnv_parser.add_argument('--sample', required=False, default='', help='Optionally, users can choose one sample to run.')
cnv_parser.set_defaults(func=discover)

args = parser.parse_args()
args.func(args)
