from hmm.Interval import *
import os

class DataLoader(object) :
    def __init__(self, filename) :
        self._filename = None 
        if os.path.exists(filename) :
            self._filename = filename
            self._datafile = open(self._filename)
            self._hasSkipHeadline = False
        
    def getTargetsList(self) :
        if self._filename :
            datafile = open(self._filename)
            headline = datafile.readline().strip('\n')
            targets = []
            #since there could be more than one chrs in a target line
            #the targets of different chrs should be put in different list
            targets_list = []
            headline_temp = headline.split('\t')[1:]
            pre_str = headline_temp[0]

            for targetstr in headline_temp:
                if targetstr.split(':')[0] == pre_str.split(':')[0]:
                    targets.append(self.buildIntervalFromStr(targetstr))
                else:
                    targets_list.append(targets)
                    targets = [self.buildIntervalFromStr(targetstr)]
                pre_str = targetstr
            targets_list.append(targets)

            return targets_list

    def buildIntervalFromStr(self, s):
        temp1 = s.split(':')
        _chr = temp1[0]
        temp2 = temp1[1].split('-')
        bp1 = temp2[0]
        bp2 = temp2[1]
        interval = Interval(_chr, bp1, bp2)
        
        return interval
        
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

    def outputCNV(self, sample, targets, pathlist, observations) :
        for i in range(len(targets)) :
		    sample.write(targets[i].getInfo() + '\t' + pathlist[i] + '\t' + observations[i] + '\n')
        
        
if __name__ == '__main__' :
    dataLoader = DataLoader("DATA.PCA_normalized.filtered.sample_zscores.RD.txt")
    print len(dataLoader.getTargets())

