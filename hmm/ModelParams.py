from numpy import *
from Interval import *
import math
from scipy import stats
import Util
from PreciseNonNegativeReal import *
import pdb

class ModelParams(object):
    state = {'DEL':0, 'DIPLOID':1, 'DUP':2}
    statestr = ['DEL', 'DIPLOID', 'DUP'] 

    def __init__(self, params, intervals, het_nums):
        self.het_nums = het_nums
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
        self._mu = []
        self._fi = []
        for i in range(3, 9, 2) :
            m = float(params[i])
            sd = float(params[i+1])
            self._mean.append(m)
            self._sd.append(sd)

        # pdb.set_trace()

        # for i in range(self.getNumHiddenStates()):
        #     self._mu.append((i+1)/ float(2) * params[-2])
        #     self._fi.append(float(params[-1]))

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
            # pdb.set_trace()
            emissProbs.append(PreciseNonNegativeReal(stats.norm.pdf(value, self._mean[i], self._sd[i])))
            # emissProbs.append(PreciseNonNegativeReal(self.calculate_NB_Prob(value,self._mu[i],self._fi[i])))
        return emissProbs

    def getEmissProbsWithSNP(self, value, t) :
        if isinstance(value, str):
            value = float(value)
        emissProbs = []

        # if self.getCH(1, t) != 1 or self.getCH(2, t) != 1:
        #     pdb.set_trace()

        for i in range(self.getNumHiddenStates()) :
            if i == 0:
                emissProbs.append(PreciseNonNegativeReal(self.getCH(0, t) * stats.norm.pdf(value, self._mean[i], self._sd[i])))
                # emissProbs.append(PreciseNonNegativeReal(self.getCH(0, t) * self.calculate_NB_Prob(value,self._mu[i],self._fi[i])))
            else:
                emissProbs.append(PreciseNonNegativeReal(self.getCH(1, t) * stats.norm.pdf(value, self._mean[i], self._sd[i])))
                # emissProbs.append(PreciseNonNegativeReal(self.getCH(1, t) * self.calculate_NB_Prob(value,self._mu[i],self._fi[i])))

        return emissProbs

    # use these probabilities to penalize the heterozgosity at the regions of deletions
    # Jixuan's original
    # def getCH(self, i, t):
    #     if i ==0 or i == 1:
    #         return 10 ** (i - 3 - self.het_nums[t])
    #     elif i >= 2:
    #         return (1 - self.getCH(0, t) - self.getCH(1, t)) / 10

    # microtan modified
    def getCH(self, i, t):
        if self.het_nums[t] > 0:
            if i == 0 :
                return 10 ** (-3 - self.het_nums[t])
            elif i >= 1:
                return (1 - self.getCH(0, t)) / 10
        else:
            return 1
    
    def getInitProbs(self):
        v = self._transitionMat[self.state['DIPLOID']]
        return v
        norm = linalg.norm(v)
        if norm==0: 
            return v
        return v/norm

    def getDiploidTransProb(self):
        return self._transitionMat[self.state['DIPLOID'], :]

    def calculate_NB_Prob(self,x,mu,fi):
        prob_nb = (math.e)**(math.lgamma(x+1/fi)-math.lgamma(x+1)-math.lgamma(1/fi)+x*math.log(mu/(mu+1/fi))-1/fi*math.log(1+mu*fi))
        return prob_nb

    @staticmethod
    def getNumHiddenStates():
        return len(ModelParams.state)

