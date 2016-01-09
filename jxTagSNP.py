from __future__ import division
import pysam
import vcf
import os,sys
import gzip
import pdb



def getTargetsByChr(filename):
    f = open(filename,'r')

    chr_targets ={}
    for line in f:
        line = line.strip()
        temp = line.split('\t')

        chrom = chrStr2Int(temp[0])
        start = int(temp[1])
        stop = int(temp[2])

        if chr_targets.has_key(chrom):
            chr_targets[chrom].append([chrom,start,stop])
        else:
            chr_targets[chrom] = [[chrom,start,stop]]

    return chr_targets


def findExonWithTagSNP(cnp_list, chr_targets, overlap_threshold):
    overlap_exon = []
    overlap_exon_unique = []

    for cnp in cnp_list:
        cnp_name = cnp[0]
        cnp_chr = chrStr2Int(cnp[1])
        cnp_start = cnp[2]
        cnp_stop = cnp[3]
        snp_pos = cnp[4]
        R2 = cnp[5]
        population = cnp[6]

        try:
            targets = chr_targets[cnp_chr]
        except KeyError:
            print 'Error: No chromosome %s in exome_region_dict' %cnp_chr

        for exon in targets:
            exon_chr = exon[0]
            exon_start = int(exon[1]) + 1 #for bed file
            exon_stop = int(exon[2]) + 1 #for bed file
            exon_size = exon_stop - exon_start + 1
            
            if exon_start > cnp_stop:
                break
            if exon_stop >= cnp_start:
                o_start = max(exon_start, cnp_start)
                o_stop = min(exon_stop, cnp_stop)

                overlap_ratio = float(o_stop - o_start + 1) / float(exon_size)

                if overlap_ratio >= overlap_threshold:
                    temp = cnp[:]
                    temp.extend(exon)
                    overlap_exon.append(temp)
                    print temp
                    #store the unique data
                    if exon not in overlap_exon_unique:
                        overlap_exon_unique.append(exon)
    
    return [overlap_exon, overlap_exon_unique]

def min(a, b):
    return a if a <= b else b

def max(a, b):
    return a if a >= b else b

def chrStr2Int(chr):
    chr = chr.replace('chr','')
    if chr == 'X' or chr == 'x':
        return 23
    elif chr == 'Y' or chr == 'y':
        return 24
    else:
        return int(chr)

def loadTagSNP(filename):
    tagsnp_f=open(filename,'r')
    cnp_dict = {}

    next(tagsnp_f) # Skip first row
    for row in tagsnp_f:
        row=row.strip()
        element = row.split('\t')
        cnp_name = element[0]
        chr = element[1].replace("chr","")
        cnp_start = int(element[2])
        cnp_stop = int(element[3])
        snp_pos = int(element[4])
        R2 = float(element[5])
        population = element[6]

        if cnp_dict.has_key(population):
            cnp_dict[population].append([cnp_name,chr,cnp_start,cnp_stop,snp_pos,R2,population])
        else:
            cnp_dict[population] = [[cnp_name,chr,cnp_start,cnp_stop,snp_pos,R2,population]]

    return cnp_dict

def findTagSNPForSample(sample_p, sample_vcf, cnp_dict):
    vcf_reader = vcf.Reader(filename=sample_vcf)
    cnp_list = []

    for cnp in cnp_dict[sample_p]:
        for vcf_snp in vcf_reader.fetch(cnp[1], cnp[4]-1, cnp[4]):  # TODO: INDEX begins from 0?
            if len(str(vcf_snp)) !=0:
        	    cnp_list.append(cnp)

    return cnp_list

