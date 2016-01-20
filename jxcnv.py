from __future__ import division
import os
import argparse
import glob
import jxcnv_functions as jf
from DataManager import *
from hmm.Model import *
from hmm.ModelParams import *
import operator 
import numpy as np
from VCFReader import *
from ParameterEstimation import *
import fileinput
import pdb

def bamlist2RPKM(args):
    try:
        import pysam
    except:
        print 'Cannot load pysam module!'
        sys.exit(0)
    try:
        # read target
        target_fn = str(args.target)
        targets = jf.loadTargets(target_fn)
        num_target = len(targets)
    except IOError:
        print 'Cannot read target file: ', target_fn
        sys.exit(0)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    try:
        bamlist_f = open(args.input)
    except IOError:
        sys.exit(0)

    for line in bamlist_f.readlines():
        line = line.strip('\n')
        temp = line.split('\t')
        sample_name = temp[0] + '.' + temp[2]
        bam_file = line.split('\t')[1]

        print 'Counting total number of reads in bam file: ', bam_file
        total_reads = float(pysam.view("-c", bam_file)[0].strip('\n'))
        print 'Found %d reads in bam file: ' % total_reads, bam_file

        f = pysam.AlignmentFile(bam_file, 'rb')

        if not f.has_index():
            print 'No index found for ', bam_file
            sys.exit(0)
    
        readcount = np.zeros(num_target)
        exon_bp = np.zeros(num_target)
        targetIDs = np.zeros(num_target)
        counter = 0

        # detect contig naming scheme here # TODO, add an optional "contigs.txt" file or automatically handle contig naming
        bam_contigs = f.references
        targets_contigs = [str(t) for t in set(map(operator.itemgetter("chr"), targets))]
                    
        targets2contigmap = {}
                    
        for targets_contig in targets_contigs:
            if targets_contig in bam_contigs:
                targets2contigmap[targets_contig] = targets_contig
            elif jf.chrInt2Str(targets_contig) in bam_contigs:
                targets2contigmap[targets_contig] = jf.chrInt2Str(targets_contig)
            elif jf.chrInt2Str(targets_contig).replace("chr","") in bam_contigs:
                targets2contigmap[targets_contig] = jf.chrInt2Str(targets_contig).replace("chr","")
            else:
                print "[ERROR] Could not find contig '%s' from %s in bam file! \n[ERROR] Perhaps the contig names for the targets are incompatible with the bam file ('chr1' vs. '1'), or unsupported contig naming is used?" % (targets_contig, target_fn)
                sys.exit(0)

        print 'Calculating RPKM values...'
        
        for t in targets:
            t_chr = targets2contigmap[str(t['chr'])]

            t_start = t['start']
            t_stop = t['stop']
            
            try:
                iter = f.fetch(t_chr, t_start, t_stop)
            except:
                print "[ERROR] Could not retrieve mappings for region %s:%d-%d. Check that contigs are named correctly and the bam file is properly indexed" % (t_chr,t_start,t_stop)
                sys.exit(0)

            for i in iter:
                if i.pos+1 >= t_start:
                    readcount[counter] += 1

            exon_bp[counter] = t_stop - t_start
            targetIDs[counter] = counter + 1
            counter += 1

        # calculate RPKM values for all targets 
        rpkm = (10**9*(readcount)/(exon_bp))/(total_reads)
        
        #rpkm_f = open(args.output+'/'+os.path.splitext(bam_file)[0].split('/')[-1]+'.rpkm', 'w')
        rpkm_f = open(args.output+'/'+sample_name+'.rpkm', 'w')
        
        rpkm_f.write('chr\tstart\tstop\trpkm\n')
        for i in range(len(rpkm)):
            rpkm_f.write(targets[i]['chr'] + '\t' + str(targets[i]['start']) + '\t' + str(targets[i]['stop']) + '\t' + str(rpkm[i]) + '\n')
        rpkm_f.close()

    bamlist_f.close()
    
def filter_rpkm(args):
    rpkm_matrix = str(args.rpkm_matrix)
    raw_dir = os.path.dirname(rpkm_matrix)
    if raw_dir != '':
        raw_dir = raw_dir + '/'

    output = rpkm_matrix
    if args.output:
        output = raw_dir + str(args.output)

    print 'Loading targets...'
    temp = jf.loadTargetsFromFirstCol(rpkm_matrix)
    targets = temp['targets']
    targets_str = temp['targets_str']

    print 'Loading parameters from ' + args.filter_params
    f = open(args.filter_params)
    params = f.readline().split('\t')
    min_gc = int(params[0])
    max_gc = int(params[1])
    min_map = int(params[2])
    max_exon = int(params[3])
    min_rpkm = float(params[4])

    #Fileter targets which median of RPKM < min_rpkm

    print "Calculating GC content..."
    #GC_percentage = jf.calGCPercentage(targets, args.ref_file)
    #jf.saveNormValues(raw_dir + 'GC_percentage', targets_str, GC_percentage, 'GC content')
    GC_percentage = jf.loadNormValues(raw_dir + 'GC_percentage')

    print 'Calculating mapping ability...'
    #map_ability = jf.calMapAbility(targets, args.map_file)
    #jf.saveNormValues(raw_dir + 'mapping_ability', targets_str, map_ability, 'Mapping ability')
    map_ability = jf.loadNormValues(raw_dir + 'map_ability')

    print 'Calculating exon length...'
    exon_length = jf.calExonLength(targets)
    #jf.saveNormValues(raw_dir + 'exon_length', targets_str, exon_length, 'Exon length')

    print 'Filtering targets by GC content, mapping ability and exon length'
    matrix = open(rpkm_matrix)
    n_lines = []
    e_lines = []
    lines = matrix.readlines()
    t = lines[0].index('\t')
    title = lines[0][0:t] + '\t' + '\t'.join(['GC Content', 'Mapping Ability', 'Exon Length']) + lines[0][t:]
    for n in range(len(lines)-1):
        # pdb.set_trace()
        line = lines[n+1]
        t = line.index('\t')
        gc = GC_percentage[n]
        m = map_ability[n]
        el = exon_length[n]
        rpkm_line = line.strip().split('\t')[1:]
        rpkm_line_median = np.median([float(s) for s in rpkm_line])
        l = line[0:t] + '\t' + '\t'.join([str(gc), str(m), str(el)]) + line[t:]
        if gc < min_gc or gc > max_gc or m < min_map or el > max_exon or rpkm_line_median < min_rpkm:
            e_lines.append(l)
        else:
            n_lines.append(l)

    n_file = open(output + '.normal', 'w')
    n_file.write(title)
    n_file.writelines(n_lines)
    n_file.close()

    e_file = open(output + '.filtered', 'w')
    e_file.write(title)
    e_file.writelines(e_lines)
    e_file.close()
                
def normalize(args): 
    rpkm_matrix = str(args.rpkm_matrix)
    raw_dir = os.path.dirname(rpkm_matrix)
    if raw_dir != '':
        raw_dir = raw_dir + '/'

    output = rpkm_matrix + '.normalized'
    if args.output:
        output = raw_dir + str(args.output) + '.normalized'

    print 'Loading matrix...'
    result = jf.loadRPKMMatrix(rpkm_matrix)
    rpkm = result['rpkm']
    samples = result['samples']
    annotation = result['annotation']
    targets = result['targets']
    targets_str = result['targets_str']

    GC_percentage = annotation[:, 0]
    map_ability = annotation[:, 1]
    exon_length = annotation[:, 2]

    GC_index = {}
    for ind in range(len(GC_percentage)):
        gc = GC_percentage[ind]

        if GC_index.has_key(gc):
            GC_index[gc].append(ind)
        else:
            GC_index[gc] = [ind]

    print 'Normalizing by GC percentage...'
    corrected_rpkm = np.zeros([len(rpkm), len(rpkm[0])], dtype=np.float)
    for i in range(len(samples)):
        print 'Normalizing RPKM by GC content for sample: ' + samples[i]
        overall_median = np.median(rpkm[:, i])
        for gc in GC_index.keys():
            t_ind = GC_index[gc]
            t_median = np.mean(rpkm[t_ind, i]) 
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, GC: %d' %(samples[i], gc)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = rpkm[t_ind, i] * overall_median / t_median

    map_index = {}
    for ind in range(len(map_ability)):
        _map = map_ability[ind]

        if map_index.has_key(_map):
            map_index[_map].append(ind)
        else:
            map_index[_map] = [ind]

    print 'Normalizing by Mapping ability...'
    for i in range(len(samples)):
        print 'Normalizing RPKM by mapping ability for sample %s' %samples[i]
        overall_median = np.median(corrected_rpkm[:, i])
        for _map in map_index.keys():
            t_ind = map_index[_map]
            t_median = np.mean(corrected_rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, Mapping ability: %d' %(samples[i], _map)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = corrected_rpkm[t_ind, i] * overall_median / t_median

    length_index = {}
    for ind in range(len(exon_length)):
        _length = exon_length[ind]
        if length_index.has_key(_length):
            length_index[_length].append(ind)
        else:
            length_index[_length] = [ind]

    print 'Normalizing by exon length...'
    for i in range(len(samples)):
        print 'Normalizing RPKM for by exon length for sample %s' %samples[i]
        overall_median = np.median(corrected_rpkm[:, i])
        for _length in length_index.keys():
            t_ind = length_index[_length]
            t_median = np.mean(corrected_rpkm[t_ind, i])
            if t_median == 0:
                print 'WARNING. Median == 0, sample: %s, Exome length: %d' %(samples[i], _length)
                corrected_rpkm[t_ind, i] = 0
            else:
                corrected_rpkm[t_ind, i] = corrected_rpkm[t_ind, i] * overall_median / t_median

    print 'Saving matrix..'
    jf.saveRPKMMatrix(output, samples, targets_str, np.transpose(corrected_rpkm))
    
def RPKM2Matrix(args):
    rpkm_dir = str(args.rpkm_dir)
    rpkm_files = glob.glob(rpkm_dir + "/*")
    if len(rpkm_files) == 0:
        print 'Cannot find any rpkm files'
        sys.exit(0)

    try:
        # read target
        target_fn = str(args.target)
        targets = jf.loadTargetsStr(target_fn)
        num_target = len(targets)
    except IOError:
        print 'Cannot read target file: ', target_fn
        sys.exit(0)
   
    samples = {}
    for f in rpkm_files:
        s = '.'.join(f.split('/')[-1].split('.')[0:-1])
        samples[s] = f

    RPKM_matrix = np.zeros([num_target, len(samples)], dtype=np.float)

    for i,s in enumerate(samples.keys()):
        t = np.loadtxt(samples[s], dtype=np.float, delimiter="\t", skiprows=1, usecols=[3])
        RPKM_matrix[:,i] = t
        print "Successfully read RPKM for" + s
    
    output = 'RPKM_matrix.raw'
    if args.output:
        output = args.output+'.raw'
    output_f = open(output, 'w')
    output_f.write('Targets\t' + '\t'.join(samples.keys()) + '\n')
    #np.savetxt(output, RPKM_matrix, fmt='%.10f', delimiter='\t', header='\t'.join(samples.keys()), comments='')
    for i in range(len(RPKM_matrix)):
        output_f.write(targets[i] + '\t' + '\t'.join(str(r) for r in RPKM_matrix[i]) + '\n')

    output_f.close()

def svd(args):
    filename = args.rpkm_matrix 
    f_dir = os.path.dirname(filename)
    if f_dir != '':
        f_dir = f_dir + '/'

    output = filename
    if args.output:
        output = f_dir + str(args.output)
   
    # count the columns number of the data file
    f = open(filename)
    temp = f.readline().strip().split('\t')
    colsnum = len(temp)
    
    # skip 1st row and 4 columns
    print 'Loading file...'
    data = np.loadtxt(filename, dtype=np.float, delimiter='\t', skiprows=1, usecols=range(4, colsnum)) 
    # loading targets str
    targets = jf.loadTargetsStrFromFirstCol(filename)
    # names of samples
    samples = temp[4:] 

    # test for conifer
    # data = np.loadtxt(filename, dtype=np.float, delimiter='\t')
    
    ### for test
    #data = np.transpose(data) 
    
    # svd transform
    # comp_removed = args.svd
    print 'SVD...'
    # pdb.set_trace()
    U, S, Vt = np.linalg.svd(data, full_matrices=False)
    
    # new_S = np.diag(np.hstack([np.zeros([comp_removed]), S[comp_removed:]]))
    index = S < 0.7 * np.mean(S)
    new_S = np.diag(S * index)
    
    # reconstruct data matrix
    data = np.dot(U, np.dot(new_S, Vt))

    # save to files
    file_u = open(output + '.U', 'w')
    file_s = open(output + '.S', 'w')
    file_vt = open(output + '.Vt', 'w')
    #file_svd = open(f_dir + args.output + '.SVD', 'w')
    
    print 'Saving SVD files...'
    np.savetxt(file_u, U, delimiter='\t')
    np.savetxt(file_s, S, delimiter='\t')
    np.savetxt(file_vt, Vt, delimiter='\t')
    #np.savetxt(file_svd, np.transpose(data), delimiter='\t')
    file_u.close()
    file_s.close()
    file_vt.close()

    print 'Saving matrix..'
    jf.saveRPKMMatrix(output + '.SVD', samples, targets, np.transpose(data))

def discover(args) :
    datafile = args.datafile
    outputfile = args.output
    paramsfile = args.params
    sample_req = args.sample
    hetsnp = args.hetsnp
    tagsnp = args.tagsnp
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

    if (hetsnp or tagsnp) and vcf_file == ''
        print 'Error: please indicate a vcf file!'
        system.exit(0)

    if vcf_file != '':
        vcf_reader = VCFReader(vcf_file)

    if tagsnp:
        cnp_dict = vcf_reader.loadTagSNP()

    while sample :
        if sample_req == '' or (sample_req != '' and sample['sample_id'] == sample_req):
            sample_flag = True

            if hetsnp:
                print 'Parsing SNV information from VCF file for: ' + sample['sample_id']
                snp_info = vcf_reader.getSNPInfo(sample['sample_id'], targets_list)
                # print snp_info
            if tagsnp:
                chr_targets = dataloader.getChrTargets()
                cnp_list = vcf_reader.findTagSNPForSample(sample['sample_pop'], cnp_dict)
                tagsnp_info = vcf_reader.findExonWithTagSNP(cnp_list, targets_list, overlap_threshold)

            #target_index is used to split observations sequence
            target_index_begin = 0
            target_index_end = 0
            temp = 1

            #estimate NB paramters from sample['observations']
            
            sample_observations = []
            remove_list = []
            if mode == 'baseline' or mode == 'reference' or mode == 'ref':
                baseline_list = ndarray.tolist(baseline_np)
                sample['observations'] = [ float(x) for x in sample['observations']]
                try:
                    sample_observations = [(v/baseline_list[i]) for i, v in enumerate(sample['observations'])]
                except Exception, e:
                    print e
                    pdb.set_trace()
                    pass
            else:
                # sample_observations  = [np.round(float(x)*10000,decimals=0) for x in sample['observations'][:-1]]
                # sample_observations  = [float(x) for x in sample['observations'][:-1]]
            
                for i,value in enumerate(sample['observations']):
                    if np.isnan(float(value)) :
                        remove_list.append(i)
                    else:
                        sample_observations.append(float(value))

            #Parameters estimation                
            # pdb.set_trace()
            # parameterLoader = ParameterEstimation(sample_observations)
            # parameterList = parameterLoader.fit(sample_observations)
            # print "Estimated Paramters: ",parameterList
            # ## print 'The parameters for negative bionomial model',parameterList
            # params_NB_mu = parameterList[0] #6.532982
            # params_NB_fi = parameterList[1] #1/3.466830 
            # params.append(params_NB_mu)
            # params.append(params_NB_fi)
            
            # if mode == 'svd' or mode == 'SVD':
            sample_observations = ndarray.tolist(stats.zscore(sample_observations))


            for targets in targets_list:
                #if targets[0]._chr == 'X' or targets[0]._chr == 'Y':
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

#BAM List -> RPKM
svd_parser = subparsers.add_parser('rpkm', help="Create RPKM matrix from a BAM list")
svd_parser.add_argument('--target', required=True, help='Target definition file')
svd_parser.add_argument('--input', required=True, help='BAM file list, each line for each sample')
svd_parser.add_argument('--output', required=True, help='Directory for RPKM files')
svd_parser.set_defaults(func=bamlist2RPKM)

#RPKM files -> Matrix
svd_parser = subparsers.add_parser('merge_rpkm', help="Merge RPKM files to a matrix")
svd_parser.add_argument('--rpkm_dir', required=True, help='RPKM files')
svd_parser.add_argument('--target', required=True, help='Target definition file')
svd_parser.add_argument('--output', required=False, help='Matrix file')
svd_parser.set_defaults(func=RPKM2Matrix)

# Filter matrix by GC content, mapping ability and exon length
svd_parser = subparsers.add_parser('filter', help="Filter matrix by GC content, mapping ability and exon length")
svd_parser.add_argument('--rpkm_matrix', required=True, help='Matrix of RPKM values')
svd_parser.add_argument('--ref_file', required=True, help='Reference file for the calculation of GC percentage')
svd_parser.add_argument('--map_file', required=True, help='Mapping ability file.')
svd_parser.add_argument('--filter_params', required=True, help='Parameters of filtering')
svd_parser.add_argument('--output', required=False, help='Filtered matrix')
svd_parser.set_defaults(func=filter_rpkm)

#normalize RPKM
svd_parser = subparsers.add_parser('norm_rpkm', help="Normalize RPKM values")
svd_parser.add_argument('--rpkm_matrix', required=True, help='Matrix of RPKM values')
svd_parser.add_argument('--output', required=False, help='Normalized RPKM matrix')
svd_parser.set_defaults(func=normalize)


#SVD
svd_parser = subparsers.add_parser('svd', help="SVD")
svd_parser.add_argument('--rpkm_matrix', required=True, help='')
svd_parser.add_argument('--output', required=False, help='')
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
cnv_parser.add_argument('--hetsnp', type=bool, required=False, default=False, help='')
cnv_parser.add_argument('--tagsnp', type=bool, required=False, default=False, help='')
cnv_parser.set_defaults(func=discover)

args = parser.parse_args()
args.func(args)
