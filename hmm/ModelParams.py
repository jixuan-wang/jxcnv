from numpy import *
from Interval import *
import math
from scipy import stats
import Util
from PreciseNonNegativeReal import *

class ModelParams(object):
    state = {'DEL':0, 'DIPLOID':1, 'DUP':2}
    statestr = ['DEL', 'DIPLOID', 'DUP'] 

    def __init__(self, params, intervals):
        self._p = float(params[0])
        self._e = int(params[1])
        self._d = int(params[2]) * Interval.KB

        self._q = 1 / float(self._e)
        self._lambda = 1 / float(self._d) 
        
        num_states = self.getNumHiddenStates()
        self._transitionMat = zeros((num_states, num_states), dtype=float)
        DEL = self.state['DEL']
        DIPLOID = self.state['DIPLOID']
        DUP = self.state['DUP']

        self._transitionMat[DEL, DEL] = 1 - self._q
        self._transitionMat[DEL, DIPLOID] = self._q
        self._transitionMat[DEL, DUP] = 0

        self._transitionMat[DIPLOID, DEL] = self._p
        self._transitionMat[DIPLOID, DIPLOID] = 1 - 2 * self._p
        self._transitionMat[DIPLOID, DUP] = self._p

        self._transitionMat[DUP, DEL] = 0 
        self._transitionMat[DUP, DIPLOID] = self._q 
        self._transitionMat[DUP, DUP] = 1 - self._q 
        
        self._mean = []
        self._sd = []
        for i in range(3, 9, 2) :
            m = float(params[i])
            sd = float(params[i+1])
            self._mean.append(m)
            self._sd.append(sd)
        
        #calculate the distances between targets which will be used to
        #calculate transition probalities
        self.interTargetDistances = []
        for i in range(len(intervals) - 1):
            self.interTargetDistances.append(intervals[i].distance(intervals[i+1]))

        return None    

    def transitionProb(self, i, j, t1):
        distInBases = self.interTargetDistances[t1] 

        f = 0
        if distInBases < float("inf") :
            f = math.exp(-(distInBases * self._lambda))

        return f * self._transitionMat[i, j] + (1-f) * self.getDiploidTransProb()[j]

    #return the transition matrix at t1->t1+1
    def transitionMatirx(self, t1) :
        tm = Util.getMatrix(self.getNumHiddenStates(), self.getNumHiddenStates())

        for i in range(self.getNumHiddenStates()) :
            for j in range(self.getNumHiddenStates()):
                tm[i][j] = PreciseNonNegativeReal(self.transitionProb(i, j, t1))

        return tm
    
    #value: read depth
    def getEmissProbs(self, value) :
        if isinstance(value, str):
            value = float(value)
        emissProbs = []
        for i in range(self.getNumHiddenStates()) :
            emissProbs.append(PreciseNonNegativeReal(stats.norm.pdf(value, self._mean[i], self._sd[i])))
        return emissProbs
        
    
    def getInitProbs(self):
        v = self._transitionMat[self.state['DIPLOID']]
        return v
        norm = linalg.norm(v)
        if norm==0: 
            return v
        return v/norm

    def getDiploidTransProb(self):
        return self._transitionMat[self.state['DIPLOID'], :]

    @staticmethod
    def getNumHiddenStates():
        return len(ModelParams.state)

