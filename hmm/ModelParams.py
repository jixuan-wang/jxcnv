class ModelParams(object):
    state = ['DEL', 'DIP', 'DUP']
    
    def __init__(self, init_probs):
        self.init_probs = init_probs
        return None    
    
    def getInitProbs(self):
        return self.init_probs


    @staticmethod
    def getNumHiddenStates():
        return len(ModelParams.state)


if __name__ == '__main__':
    mp = ModelParams([0.2, 0.7, 0.1])
    print mp.getInitProbs()[0]
