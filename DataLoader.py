from hmm.Interval import *
import os

class DataLoader(object) :
    def __init__(self, filename) :
        self._filename = None 
        if os.path.exists(filename) :
            self._filename = filename
            self._datafile = open(self._filename)
            self._hasSkipHeadline = False
        
    def getTargets(self) :
        if self._filename :
            datafile = open(self._filename)
            headline = datafile.readline().strip('\n')
            targets = []
            for targetstr in headline.split('\t')[1:] :
                temp1 = targetstr.split(':')
                _chr = temp1[0]
                temp2 = temp1[1].split('-')
                bp1 = temp2[0]
                bp2 = temp2[1]
                interval = Interval(_chr, bp1, bp2)
                targets.append(interval)
            return targets

    def skipHeadline(self):
        if self._datafile and not self._hasSkipHeadline :
            self._datafile.readline()
            self._hasSkipHeadline = True

    def getNextSample(self) :
        samplestr = self._datafile.readline()
        if samplestr :
            samplestr.strip('\n')
            samlist = samplestr.split('\t')
            sample = {'sample_id':samlist[0], 'observations':samlist[1:]}

            return sample

        return None

    def getParams(self, paramsfile):
        try:
            pfile = open(paramsfile)
            pline = pfile.readline().strip('\n')
            params = pline.split('\t')
            pfile.close()
            
            return params
        except IOError, msg:
            print msg
            return None
        
    def closeFile(self):
        if self._datafile:
            self,_datafile.close()

    def outputCNV(self, sample, targets, pathlist) :
        output = file(sample, 'w')

        for i in range(len(targets)) :
		    output.write(targets[i].getInfo() + pathlist[i] + '\n')

        output.close()
        
        
if __name__ == '__main__' :
    dataLoader = DataLoader("DATA.PCA_normalized.filtered.sample_zscores.RD.txt")
    print len(dataLoader.getTargets())

