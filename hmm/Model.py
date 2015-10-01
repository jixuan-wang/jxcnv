from numpy import *

class Model (object)
    def __init__(self, _model_params, _observations):
	    self._model_params = _model_params
        self._observations = _observations

    def forwardBackward_Viterbi(self):
        s = _model_params.getNumHiddenStates()
        o = len(self._observations)
        alpha = mat(zeros((o, s)))
        beta = mat(zeros((o, s)))
        gamma = mat(zeros((o, s)))
        
        # compute forward probabilities
        for t in range(o):
             print 'Calculating forward message ' + str(t)
             calcForwardMessage(t, alpha)
             print alpha[t]
        
        #compute backward probabilities
        for t in range(o)[::-1] :
            print 'Calculating backward message ' + str(t)
            calcBackwardMessage(t, beta)
            print beta[t]
        
        for t in range(o) :
            gamma[t] = multiply(alpha[t], beta[t])

        _vpath = []
        for t in range(o) :
            _vpath.append(argmax(gamma[t]))

        return _vpath

    def calcForwardMessage(self, t, alpha):
        forwardMess = alpha[t]
        s = self._model_params.getNumHiddenStates()

        if t == 0:
            for i in range(s):
                forwardMess[i] = self._model_params.getInitProbs()[i]
        else:
            incomingForwardMess = alpha[t-1]
            transMat = self._model_params.transitionMatirx(t-1)

            forwardMess = incomingForwardMess * transMat
        
        multiplyMessageTimesEmissProbs(forwardMess, t) 

    def calcBackwardMessage(self, t, beta) :
        backwardMess = beta[t]
        s = self._model_params.getNumHiddenStates()

        if t == s - 1 :
            backwardMess = mat(ones((1, s)))  
        else :
            incomingBackwardMess = bete[t+1]
            multiplyMessageTimesEmissProbs(incomingBackwardMess, t+1)
            transMatTranspose = transpose(self._model_params.transitionMatirx(t))

            backwardMess = incomingBackwardMess * transMatTranspose

            
    def multiplyMessageTimesEmissProbs(self, message, t):
        emissiProbs = self._model_params.getEmissiProbs(self._observation[t])

        for i in range(self._model_params.getNumHiddenStates()) :
            message[0, i] *= emissProbs[i]
        








        


