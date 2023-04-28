import numpy as np
from scipy.stats import wasserstein_distance as WashMyBalls

x1 = np.random.rand(1,10)
x2 = np.random.rand(1,10)
x3 = np.random.rand(1,10)
x4 = np.random.rand(1,10)
x5 = np.random.rand(1,10)
x6 = np.random.rand(1,10)
x7 = np.random.rand(1,10)
x8 = np.random.rand(1,10)
x9 = np.random.rand(1,10)
x10 = np.random.rand(1,10)

SampleData = [x1,x2,x3,x4,x5,x6,x7,x8,x9,x10]

def WasserSwag(Data):
    IndexPairs = []
    Total = 0
    for i in range(len(Data)-1):
        for n in range(i+1, len(Data)):
            if all(elem in IndexPairs for elem in [[i, n], [n, i]]):
                continue
            IndexPairs.append([i,n])
            D1 = Data[i][0].tolist()
            D2 = Data[n][0].tolist()
            WashYourBalls = WashMyBalls(D1,D2)
            Total += WashYourBalls
    return Total
    
print(WasserSwag(SampleData))






    
