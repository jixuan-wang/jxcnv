from PreciseNonNegativeReal import *

def getMatrix(x, y, default=PreciseNonNegativeReal(0):
    m = []
    for i in range(0, 3):
        tmp = []
        for j in range(0, 3):
            tmp.append(default)
        m.append(tmp)
