import numpy as np
import pysam
import pdb
import bx.bbi.bigwig_file

def loadMatrixFromFile(filename, skiprows, skipcols, type='float', delimiter='\t'):
    return

def loadTargets(target_filename):
    targetfile = open(target_filename)
    targets = []
    target_id = 1

    for line in targetfile.readlines():
        line = line.strip('\n')
        #temp1 = line.split(':')
        #temp2 = temp1[1].split('-')
        #targets.append({'targetID': target_id, 'chr': temp1[0], 'start': temp2[0], 'stop': temp2[1]})
        #target_id += 1
        temp = line.split('\t')
        targets.append({'targetID': int(temp[3]), 'chr': temp[0], 'start': int(temp[1]), 'stop': int(temp[2])})


    return targets

def loadTargetsStr(target_filename):
    targetfile = open(target_filename)
    targets = []

    for line in targetfile.readlines():
        line = line.strip()
        temp = line.split('\t')
        targets.append(temp[0]+':'+temp[1]+'-'+temp[2])

    return targets

def loadTargetsFromFirstCol(filename):
    f = open(filename)
    lines = f.readlines()
    targets = []
    for line in lines[1:len(lines)]:
        chr = line.split('\t')[0].split(':')[0]
        start = int(line.split('\t')[0].split(':')[1].split('-')[0])
        stop = int(line.split('\t')[0].split(':')[1].split('-')[1])
        targets.append({'chr': chr, 'start': start, 'stop': stop})
    return targets

def loadRPKMMatrix(filename):
    f = open(filename)
    line = f.readline()
    line = line.strip()
    colsnum = len(line.split('\t'))

    data = np.loadtxt(f, dtype=np.float, delimiter='\t', skiprows=0, usecols=range(1, colsnum)) 
    return {'samples': line.split('\t')[1:], 'data': data}
    

def chrInt2Str(chromosome_int):
    return 'chr' + str(chromosome_int)
    if int(chromosome_int) == 23:
        return 'chrX'
    elif int(chromosome_int) == 24:
        return 'chrY'
    else:
        return 'chr' + str(chromosome_int)

def calGCPercentage(targets, ref_file):
    fasta_f = pysam.FastaFile(ref_file)
    GC_percentage = []

    for i in range(len(targets)):
        r_region = fasta_f.fetch(targets[i]['chr'], targets[i]['start'], targets[i]['stop'])
        reg_len = targets[i]['stop'] - targets[i]['start'] + 1
        GC = 0
        n_num = 0

        for b in range(len(r_region)):
            if str(r_region[b]).upper() == 'G' or str(r_region[b]).upper() == 'C':
                GC += 1
            elif str(r_region[b]).upper() == 'N':
                n_num += 1

        if n_num / reg_len > 0.2:
            GC_p = -1
        else:
            GC_p = GC * 100 / reg_len

        GC_percentage.append([i, GC_p])
        print i, GC_p
    
    return GC_percentage

def calMapAbility(targets, map_file):
    bw = bx.bbi.bigwig_file.BigWigFile(open(map_file, "rb"))
    map_ability = []
    
    for i in range(len(targets)):
        t = targets[i] 
        map_summary = bw.query(chrInt2Str(t['chr']), t['start'], t['stop'], 1)
        try:
            _map = int(map_summary[0]['mean']*100+0.5)
        except:
            _map = -1

        map_ability.append([i, _map])
        print i, _map
    return map_ability

def calExonLength(targets):
    exon_length = []
    for i in range(len(targets)):
        t = targets[i]
        l = np.round((t['stop']-t['start']+1)/50)*50
        exon_length.append([i, l])
    return exon_length


def groupBy(data):
    result = {}
    exclude = []
    for ind, val in data:
        if val == -1:
            exlude.append(ind)
        elif result.has_key(val):
            result[val].apend(ind)
        else:
            result[val] = [ind]
    return {'exlude': exclude, 'dict': result}
