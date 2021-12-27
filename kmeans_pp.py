import pandas as pd
import numpy as np

"""names = []
DataPoints = pd.read_csv(r'C:\Users\weamm\Downloads\input_1_db_1.txt' , names =names + ["dim" + str(i) for i in range(1,5)])
DataPoints["D"] = [0 for i in range(len(DataPoints))]
DataPoints["P"] = [0 for i in range(len(DataPoints))]
print(DataPoints)
print(DataPoints)"""


def dist(point1, point2):
    sum = 0
    d = len(point1)
    for i in range(d):
        sum += (point1[i] - point2[i]) ** 2
    sum = sum ** 0.5
    return sum


def mindist(centroids_index, point2, DataPoints,dimension):
    min = "inf"
    index = 0

    for i in range(len(centroids_index)):
        point1 = []
        for j in range(1, dimension+1):
            point1.append(DataPoints["dim"+str(j)][centroids_index[i]])
        distance = dist(point1, point2)
        if (min == "inf" or distance < min):
            min = distance
            index = i

    return centroids_index[index]


def kmeansPlus(k, max_iter, epsilon, filename1, filename2):
    centroids = []
    centroids_index = []


def init(k, filename1, filename2, centroids, centroids_index):

    DataPoints = pd.read_csv(filename1)
    dimension = len(DataPoints.column)
    DataPoints.columns = ["dim" + str(i) for i in range(dimension)]

    DataPoints["D"] = [0 for i in range(len(DataPoints))]
    DataPoints["P"] = [0 for i in range(len(DataPoints))]

    np.random.seed(0)
    muo = np.random.choice(len(DataPoints), 1)
    for i in range(dimension):
        point1 = [DataPoints["dim" + str(i)][muo]]

    print(DataPoints)

    """for i in range len(DataPoints):
        point1 = 
        DataPoints["D"][i] = dist(DataPoints)"""

    """i=1
    while(i<k):
        muo = np.random.choice(len(DataPoints) , 1, p=DataPoints['P'].to_numpy())"""



