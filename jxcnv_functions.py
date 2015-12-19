import numpy

def loadMatrixFromFile(filename, skiprows, skipcols, type='float', delimiter='\t'):

def loadTargets(target_filename):
    targetfile = open(target_filename)
    targets = []
    target_id = 1

    while line = targetfile.readline():
        line = line.strip('\n')
        #temp1 = line.split(':')
        #temp2 = temp1[1].split('-')
        #targets.append({'targetID': target_id, 'chr': temp1[0], 'start': temp2[0], 'stop': temp2[1]})
        #target_id += 1
        temp = line.split('\t')
        targets.append({'targetID': temp[3], 'chr': temp[0], 'start': temp[1], 'stop': temp[2]})


    return targets

def chrInt2Str(chromosome_int):
    if int(chromosome_int) == 23:
        return 'chrX'
    elif int(chromosome_int) == 24:
        return 'chrY'
    else:
        return 'chr' + str(chromosome_int)

        



