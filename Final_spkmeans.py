import pandas as pd
import numpy as np
import mykmeanssp
import argparse
import sys


def dist(point1, point2):
    sum = 0
    d = len(point1)
    for i in range(d):
        sum += (point1[i] - point2[i]) ** 2

    return sum


def kmeansPlus(k, goal, filename):
    DataPoints = ReadData(k, filename)
    Data_indices = DataPoints.iloc[:, :1]
    Data_indices = Data_indices.to_numpy()
    Data_indices = np.array(Data_indices)

    # DataPoints = DataPoints.iloc[:, 1:]
    DataPoints = DataPoints.to_numpy()
    DataPoints = np.array(DataPoints)

    rows = DataPoints.shape[0]
    dimension = DataPoints.shape[1]
    if (goal != "spk"):
        if (goal == "wam"):
            goal = 2
        if (goal == "ddg"):
            goal = 3
        if (goal == "lnorm"):
            goal = 4
        if (goal == "jacobi"):
            goal = 5
        Data_list = DataPoints.tolist()
        matrix = mykmeanssp.fit(Data_list, None, rows, dimension, k, 300, goal, 0)
        matrix = np.array(matrix)
        matrix = np.round(matrix, decimals=4)
        print_points(matrix)
        print()
        return
    else:
        K=k
        goal = 1
        if (k >= rows):
            print("Invalid Input!")
            return -1

        Data_list = DataPoints.tolist()
        if(K==0):
            matrix = mykmeanssp.fit(Data_list, None, rows, dimension, k, 300, goal, 0)
            matrix = np.array(matrix)
            matrix = np.round(matrix, decimals=4)
            K = matrix[0][0]

        matrix = mykmeanssp.fit(Data_list, None, rows, dimension, K, 300, goal, 0)
        matrix = np.array(matrix)
        matrix = np.round(matrix, decimals=4)

        centroids = np.ndarray((K, K), float)
        centroids_index = np.ndarray(K, int)
        init_Centroids(matrix, centroids, centroids_index, K, K, rows)
        matrix = matrix.tolist()
        centroids_list = centroids.tolist()
        centroids = mykmeanssp.fit(matrix, centroids_list, rows, K, K, 300, goal, 1)
        centroids = np.array(centroids)
        centroids = np.round(centroids, decimals=4)
        print_points(centroids)
        print()
        return
    return


def print_points(centroids):
    for i in range(len(centroids)):
        centroid = centroids[i]
        for j in range(len(centroid)):
            if (j != (len(centroid) - 1)):
                print(str(centroid[j]) + ",", end="")
            else:
                if (i == len(centroids) - 1):
                    print(centroid[j], end="")
                else:
                    print(centroid[j])


def ReadData(k, filename):
    DataPoints = pd.read_csv(filename, header=None)
    return DataPoints


def init_Centroids(DataPoints, centroids, centroids_index, k, dimension, rows):
    sum1 = 0
    D = np.zeros(rows)
    P = np.zeros(rows)

    np.random.seed(0)
    mew = np.random.choice(rows)
    centroids_index[0] = mew
    centroids[0] = DataPoints[mew]

    Z = 1
    while Z < k:
        for i in range(0, rows):
            min = float("inf")
            for j in range(0, Z):
                distance = dist(DataPoints[i], centroids[j])
                if (distance < min):
                    min = distance
            sum1 -= D[i]
            D[i] = min
            sum1 += D[i]

        P = np.divide(D, sum1)
        index1 = np.random.choice(rows, p=P)
        centroids_index[Z] = index1
        centroids[Z] = DataPoints[index1]
        Z += 1
    print(','.join(str(i) for i in centroids_index), flush=True)


def start():  # gets arguments and starts the algorithm.
    parser = argparse.ArgumentParser()
    parser.add_argument("K", type=int, help="K is the number of clusters")
    parser.add_argument("goal", type=str, help="The goal which is need to calculate")
    parser.add_argument("file_name", type=str, help="The path to file which contains N observations")
    args = parser.parse_args()
    K = args.K
    goal = args.goal
    file_name = args.file_name
    # assertions

    if K == None:
        print("Invalid Input!")
        return -1
    if (K < 0):
        print("Invalid Input!")
        return -1
    goals = ["wam", "ddg", "lnorm", "spk", "jacobi"]
    if (goal not in goals):
        print("Invalid Input!")
        return -1
    if (file_name) == None:
        print("Invalid Input!")
        return -1
    kmeansPlus(K, goal, file_name)


start()
