import numpy as np

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

SampleData = np.array([x1,x2,x3,x4,x5,x6,x7,x8,x9,x10])
IndexPairs = []
for i in range(len(SampleData)-1):
    Array1 = SampleData[i]
    for n in range(i+1, len(SampleData)):
        Array2 = SampleData[n]
        if all(elem in IndexPairs for elem in [[i, n], [n, i]]):
            continue
        IndexPairs.append([i,n])
        







    
