import pandas as pd
import numpy as np


"""DataPoints = pd.read_csv(r'C:\Users\weamm\Downloads\input_1_db_1.txt')
dimension = len(DataPoints.columns)
DataPoints.columns = ["dim" + str(i) for i in range(1,dimension+1)]
DataPoints["D"] = [0 for i in range(len(DataPoints))]
DataPoints["P"] = [0 for i in range(len(DataPoints))]
point1=[]
for i in range (1,dimension+1):
        point1.append(DataPoints["dim" + str(i)][2])
print(point1)"""

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

    return min

def kmeansPlus(k,max_iter, epsilon , filename1, filename2):
    centroids = []
    centroids_index=[]



def init(k, filename1 , filename2 , centroids, centroids_index):

    sum1 = 0
    DataPoints = pd.read_csv(filename1)
    dimension = len(DataPoints.column)
    DataPoints.columns = ["dim" + str(i) for i in range(1,dimension+1)]
    
    DataPoints["D"] = [0 for i in range(len(DataPoints))]
    DataPoints["P"] = [0 for i in range(len(DataPoints))]
    
    np.random.seed(0)
    muo = np.random.choice(len(DataPoints), 1)
    
    centroids_index.append(muo)
    
    point1=[]
    for i in range (1,dimension+1):
        point1.append(DataPoints["dim" + str(i)][muo])
    
        
    for j in range(len(DataPoints)):
        point2=[]
        for i in range (1,dimension+1):
            point2.append(DataPoints["dim" + str(i)][j])
            
        newD = dist(point1 , point2)
        oldD = DataPoints["D"][j]
        DataPoints["D"][j] = newD
        sum1 = sum1 - oldD  + newD
        
    for j in range(len(DataPoints)):
        DataPoints["P"][j] = DataPoints["D"][j]/sum1
        
    i=1
    while(i<k):
        muo = np.random.choice(len(DataPoints) , 1, p=DataPoints['P'].to_numpy())
        centroids_index.append(muo)
    
        for j in range(len(DataPoints)):
            point=[]
            for i in range (1,dimension+1):
                point.append(DataPoints["dim" + str(i)][j])
            
            D = mindist(centroids_index, point, DataPoints, dimension)
            sum1 = sum1 - DataPoints["D"][j] + D
            DataPoints["D"][j] = D
            
        for j in range(len(DataPoints)):
            DataPoints["P"][j] = DataPoints["D"][j]/sum1


        
        


