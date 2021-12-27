import pandas as pd
import numpy as np



def dist(point1, point2):
    sum = 0
    d = len(point1)
    for i in range(d):
        sum += (point1[i] - point2[i]) ** 2
    sum = sum ** 0.5
    return sum


def mindist(centroids, point2, DataPoints, dimension, z):
    min = "inf"

    for i in range(z+1):
        point1 = centroids[i]
        distance = dist(point1, point2)
        if (min == "inf" or distance < min):
            min = distance

    return min


def kmeansPlus(k, MAX_ITER, filename1, filename2):
    DataPoints = ReadData(k,filename1, filename2)

    Data_indices = DataPoints.iloc[:, :1]
    Data_indices = Data_indices.to_numpy()
    Data_indices = np.array(Data_indices)

    DataPoints = DataPoints.iloc[:, 1:]
    DataPoints = DataPoints.to_numpy()
    DataPoints = np.array(DataPoints)

    rows = DataPoints.shape[0]
    dimension = DataPoints.shape[1]

    if (k >= rows):
        print("Invalid Input!")
        return -1

    centroids = np.ndarray((k, dimension), float)
    centroids_index = np.ndarray(k, int)

    init_Centroids(DataPoints, centroids, centroids_index, k, dimension, rows  )
    print(centroids_index)
    print(centroids)


def ReadData(k, filename1, filename2):
    DataPoints1 = pd.read_csv(filename1, header=None)
    DataPoints2 = pd.read_csv(filename2, header=None)
    DataPoints = pd.merge(DataPoints1, DataPoints2, on=0)
    DataPoints = DataPoints.sort_values(by=[0])

    return DataPoints


def init_Centroids(DataPoints, centroids, centroids_index, k, dimension, rows):
    sum1 = 0
    D = np.zeros(rows)
    P = np.zeros(rows)

    np.random.seed(0)
    mew = np.random.choice(rows)
    centroids_index[0] = mew
    centroids[0] = DataPoints[mew]

    point1 = DataPoints[mew]
    for j in range(rows):
        point2 = DataPoints[j]
        D[j] = dist(point1, point2)
        sum1 = sum1 + D[j]

    P = np.divide(D, sum1)

    i = 1
    while (i < k):
        mew = np.random.choice(rows, p=P)
        centroids_index[i] = mew
        centroids[i] = DataPoints[mew]

        for j in range(rows):
            point = DataPoints[j]
            d = mindist(centroids, point, DataPoints, dimension, i)

            sum1 = sum1 - D[j] + d
            D[j] = d

        P = np.divide(D, sum1)
        i+=1

import sys
args = sys.argv
maxiter = 0

if len(args) == 5:
    if (args[1].isdigit() == False or args[2].isdigit() == False):
        print("Invalid Input!")
        sys.exit()

    assert args[1].isdigit(), "Invalid input."
    assert args[2].isdigit(), "Invalid input."
    k = int(args[1])
    maxiter = int(args[2])
    inputfile1 = args[3]
    inputfile2 = args[4]
    kmeansPlus(k, maxiter, inputfile1, inputfile2)

if len(args) == 4:
    if (args[1].isdigit() == False):
        print("Invalid Input!")
        sys.exit()

    k = int(args[1])
    maxiter = 200
    inputfile = args[2]
    outputfile = args[3]
    kmeansPlus(k, maxiter, inputfile1, inputfile2)



