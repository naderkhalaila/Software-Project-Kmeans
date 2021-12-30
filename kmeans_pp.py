import pandas as pd
import numpy as np
import mykmeanssp


def dist(point1, point2):
    sum = 0
    d = len(point1)
    for i in range(d):
        sum += (point1[i] - point2[i]) ** 2
    
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
    centroids =  mykmeanssp.fit(DataPoints, centroids, rows, dimension, k, MAX_ITER)
    centroids = np.array(centroids)
    centroids = np.round(centroids, decimals = 4)
    print_centroids(centroids)
    return None


def print_centroids(centroids):
    for i in range (len(centroids)):
        centroid = centroids[i]
        for j in range(len(centroid)):
            if(j != (len(centroid)-1)):
                print(str(centroid[j])+ ",", end="")
            else:
                if(i==len(centroids)-1):
                    print(centroid[j], end="")
                else:
                    print(centroid[j])



def ReadData(k, filename1, filename2):
    DataPoints1 = pd.read_csv(filename1, header=None)
    DataPoints2 = pd.read_csv(filename2, header=None)
    DataPoints = pd.merge(DataPoints1, DataPoints2, on=0)
    DataPoints = DataPoints.sort_values(by=[0])

    return DataPoints


def init_Centroids(DataPoints, centroids, centroids_index, k, dimension, rows):
    sum1=0
    D = np.zeros(rows)
    P = np.zeros(rows)

    np.random.seed(0)
    mew = np.random.choice(rows)
    centroids_index[0] = mew
    centroids[0] = DataPoints[mew]

    Z = 1
    while Z<k:
        for i in range(0, rows):
            min = float("inf")
            for j in range(0, Z):
                distance = dist(DataPoints[i] , centroids[j])
                if(distance<min):
                    min = distance
            sum1-=D[i]
            D[i] = min
            sum1+=D[i]

        P = np.divide(D, sum1)
        index1 = np.random.choice(rows, p=P)
        centroids_index[Z] = index1
        centroids[Z] = DataPoints[index1]
        Z+=1



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



ff = kmeansPlus(7 , 300 ,r'C:\Users\weamm\Downloads\input_2_db_1.txt' , r'C:\Users\weamm\Downloads\input_2_db_2.txt' )
