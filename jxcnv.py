from __future__ import division
import argparse
from DataManager import *
from hmm.Model import *
from hmm.ModelParams import *
import numpy as np
from VCFReader import *
from ParameterEstimation import *
import pdb

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
    vcf_file = args.vcf
    mode = args.mode
    sample_flag = False #used to check whether sample_req exists

    # Build a reference set 

    if mode == 'baseline' or mode == 'reference' or mode == 'ref':
        print 'Building the reference dataset...'
        dataloader = DataManager(datafile)
        samples_np = dataloader.getAllSamples()
        baseline_np = np.average(samples_np, axis=0)
        dataloader.closeFile()
        print 'Baseline is Done.'

    print 'Loading data file...'
    dataloader = DataManager(datafile)
    params = dataloader.getParams(paramsfile)
    dataloader.skipHeadline()
    sample = dataloader.getNextSample()

    targets_list = dataloader.getTargetsList()
    output_aux = file(outputfile+'.aux', 'w')
    output_aux.write('SAMPLE_ID\tCNV_TYPE\tFULL_INTERVAL\tINDEX\tINTERVAL\tREAD_DEPTH\n')
    output = file(outputfile,'w')
    output.write('SAMPLE_ID\tCNV_TYPE\tINTERVAL\tCHROMOSOME\tSTART\tSTOP\tLENGTH\n')

    vcf_reader = None
    if vcf_file != '':
        vcf_reader = VCFReader(vcf_file)
     
    while sample :
        if sample_req == '' or (sample_req != '' and sample['sample_id'] == sample_req):
            sample_flag = True
            snp_info = list()
            if vcf_reader:
                print 'Parsing SNV information from VCF file for: ' + sample['sample_id']
                snp_info = vcf_reader.getSNPInfo(sample['sample_id'], targets_list)
                # print snp_info
            #target_index is used to split observations sequence
            target_index_begin = 0
            target_index_end = 0
            temp = 1

            #estimate NB paramters from sample['observations']
            
            sample_observations = []
            remove_list = []
            if mode == 'baseline' or mode == 'reference' or mode == 'ref':
                baseline_list = ndarray.tolist(baseline_np)
                sample['observations'][:-1] = [ float(x) for x in sample['observations'][:-1]]
                sample_observations = [(v/baseline_list[i]) for i, v in enumerate(sample['observations'][:-1])]
            else:
                # sample_observations  = [np.round(float(x)*10000,decimals=0) for x in sample['observations'][:-1]]
                # sample_observations  = [float(x) for x in sample['observations'][:-1]]
            
                for i,value in enumerate(sample['observations'][:-1]):
                    if np.isnan(float(value)) :
                        remove_list.append(i)
                    else:
                        sample_observations.append(float(value))

            #Parameters estimation                
            # parameterLoader = ParameterEstimation(sample_observations)
            # parameterList = parameterLoader.fit(sample_observations)
            # print "Estimated Paramters: ",parameterList
            ## print 'The parameters for negative bionomial model',parameterList
            ## params_NB_mu = parameterList[0] #6.532982
            ## params_NB_fi = parameterList[1] #1/3.466830 
            ## params.append(params_NB_mu)
            ## params.append(params_NB_fi)

            for targets in targets_list:
                print 'Running HMM for sample[' + sample['sample_id'] + ']: ',
                print 'chr' + targets[0]._chr + ' [' + str(temp) + '|' + str(len(targets_list)) + ']'

                temp += 1
                target_index_end = target_index_begin + len(targets)

                #Run the HMM 
                modelParams = ModelParams(params, targets, snp_info[target_index_begin:target_index_end])
                #the 'observations' of sample is splitted
                # model = Model(modelParams, sample['observations'][target_index_begin:target_index_end])
                model = Model(modelParams, sample_observations[target_index_begin:target_index_end])
                pathlist = list()
                
                if vcf_reader:
                    pathlist = model.forwardBackward_Viterbi(if_snp = True)
                else:
                    pathlist = model.forwardBackward_Viterbi(if_snp = False)
                # pdb.set_trace()
                dataloader.outputCNVaux(output_aux, sample['sample_id'], targets, pathlist, sample_observations[target_index_begin:target_index_end])
                dataloader.outputCNV(output, sample['sample_id'], targets, pathlist, sample_observations[target_index_begin:target_index_end])

                target_index_begin = target_index_end
        sample = dataloader.getNextSample()

    output.close()
    output_aux.close()
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
cnv_parser.add_argument('--mode',required=True, default='SVD', help='Data normalization by SVD or baseline mode.')
cnv_parser.add_argument('--output', required=True, help='Output file.')
cnv_parser.add_argument('--sample', required=False, default='', help='Optionally, users can choose one sample to run.')
cnv_parser.add_argument('--vcf', required=False, default='', help='Optionally, users can input snp information by specifing a vcf file')
cnv_parser.set_defaults(func=discover)

args = parser.parse_args()
args.func(args)
