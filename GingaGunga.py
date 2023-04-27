import numpy as np
from scipy.stats import wasserstein_distance as WashMyBalls
import random

x1 =[1,2,3]
x2 =[4,5,6]

""" 
x1 = list(np.random.rand(1,10))
x2 = list(np.random.rand(1,10)) """
""" x3 = list(np.random.rand(1,10))
x4 = np.random.rand(1,10)
x5 = np.random.rand(1,10)
x6 = np.random.rand(1,10)
x7 = np.random.rand(1,10)
x8 = np.random.rand(1,10)
x9 = np.random.rand(1,10)
x10 = np.random.rand(1,10) """

# SampleData = [x1,x2]

print(WashMyBalls(x1,x2))


""" IndexPairs = []
for i in range(len(SampleData)-1):
    Array1 = SampleData[i]
    for n in range(i+1, len(SampleData)):
        Array2 = SampleData[n]
        if all(elem in IndexPairs for elem in [[i, n], [n, i]]):
            continue
        IndexPairs.append([i,n])
         """


Total = WashMyBalls(x1,x2)


def WasserSwag(Data):
    IndexPairs = []
    Total = 0
    for i in range(len(Data)-1):
        for n in range(i+1, len(Data)):
            if all(elem in IndexPairs for elem in [[i, n], [n, i]]):
                continue
            IndexPairs.append([i,n])
            print(i,n)
        D1 = Data[i]
        D2 = Data[n]
        print(D1,D2)
        WashYourBalls = WashMyBalls(D1,D2)
    return Total
    






    
