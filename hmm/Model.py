from numpy import *

class Model (object)
    def __init__(self, _model_params):
	    self._model_params = _model_params

    def forwardBackward_Viterbi(self, observations):
        s = _model_params.getNumHiddenStates()
        o = len(observations)
        alpha = mat(zeros((o, s)));
        beta = mat(zeros((o, s)));
        
        # compute forward probabilities
        for t in range(o):
             print 'Calculating forward message ' + str(t)
             calcForwardMessage(t, alpha)
             print alpha[t]

    
    def calcForwardMessage(t, alpha):
        forwardMess = alpha[t]
        s = _model_params.getNumHiddenStates()

        if t == 0:
            for i in range(s):
                forwardMess[i] = _model_params.getInitProbs()[i]
        else:
            
        


