from numpy import *
from ModelParams import *
import pdb
import copy

class Model (object) :
    def __init__(self, _model_params, _observations):
        self._model_params = _model_params
        self._observations = _observations

    def forwardBackward_Viterbi(self):
        s = ModelParams.getNumHiddenStates()
        o = len(self._observations)
        alpha = mat(zeros((o, s)))
        beta = mat(zeros((o, s)))
        gamma = mat(zeros((o, s)))
        
        # compute forward probabilities
        for t in range(o):
             print 'Calculating forward message ' + str(t),
             self.calcForwardMessage(t, alpha)
             print alpha[t]
        
        #compute backward probabilities
        for t in range(o)[::-1] :
            print 'Calculating backward message ' + str(t),
            self.calcBackwardMessage(t, beta)
            print beta[t]
        
        for t in range(o) :
            gamma[t] = multiply(alpha[t], beta[t])

        _vpath = []
        for t in range(o) :
            _vpath.append(argmax(gamma[t]))

        vpathlist = []
        i = 0
        for v in _vpath :
            vpathlist.append(ModelParams.statestr[v])
            i += 1
        return vpathlist

    def calcForwardMessage(self, t, alpha):
        forwardMess = alpha[t]
        s = ModelParams.getNumHiddenStates()

        if t == 0:
            for i in range(s):
                forwardMess[0,i] = self._model_params.getInitProbs()[i]
        else:
            incomingForwardMess = alpha[t-1]
            transMat = self._model_params.transitionMatirx(t-1)

            forwardMess = incomingForwardMess * transMat
        
        self.multiplyMessageTimesEmissProbs(forwardMess, t)
        alpha[t] = forwardMess

    def calcBackwardMessage(self, t, beta) :
        backwardMess = beta[t]
        o = len(self._observations)
        s = ModelParams.getNumHiddenStates()

        if t == o - 1 :
            backwardMess = mat(ones((1, s)))  
        else :
            incomingBackwardMess = copy.deepcopy(beta[t+1])
            self.multiplyMessageTimesEmissProbs(incomingBackwardMess, t+1)
            transMatTranspose = transpose(self._model_params.transitionMatirx(t))

            backwardMess = incomingBackwardMess * transMatTranspose

        beta[t] = backwardMess
            
    def multiplyMessageTimesEmissProbs(self, message, t):
        emissProbs = self._model_params.getEmissProbs(self._observations[t])

        for i in range(ModelParams.getNumHiddenStates()) :
            message[0, i] *= emissProbs[i]

        








        


