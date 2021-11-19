def dist(point1, point2):
    
    sum=0
    d = len(point1)
    for i in range(d):
        sum+= (point1[i] - point2[i])**2
    sum = sum**0.5
    return sum

def mindist(centroids_list ,point):
    min = "inf"
    index = 0 
    for j in range (len(centroids_list)):
        distance= dist (centroids_list[j] , point)
        if (min == "inf" or distance < min):
            min = distance
            index = j 
    return index


def clustering(centroids_list ,DataPoints, indexingList):
    for i in range(len(DataPoints)):
        x = mindist(centroids_list ,DataPoints[i])
        indexingList[x].append(i)
        
      
def update(centroids_list ,DataPoints, indexingList):
    boolean = True
    for i in range(len(centroids_list)):
        old_avg = centroids_list[i]
        centroids_list[i]= Avg(indexingList[i] , DataPoints)
        distance = dist (centroids_list[i] , old_avg)
        if (distance>0.001):
            boolean = False
    
    return boolean


def Avg(lst, DataPoints):
    avg_point=[]
    d= len(DataPoints[0])
    for i in range(d):
        sum = 0 
        for item in lst:
            sum+= DataPoints[item][i]
        sum = sum/len(lst)
        avg_point.append(sum)
    return avg_point
    

def kmeans (k ,max_iter = 200, inputfile=None , outputfile=None):
    
    assert type(k) == int , "Invalid input."
    assert type(max_iter) == int , "Invalid input."
    
    centroids_list = []
    DataPoints = []

    file = open(inputfile)
    content = file.read()
    DataPoints = content.splitlines()
    for i in range(len(DataPoints)):
        DataPoints[i] = DataPoints[i].split(",")
        for j in range(len(DataPoints[i])):
            DataPoints[i][j] = float(DataPoints[i][j])
            
    for i in range(k):
        centroids_list.append(DataPoints[i])
    
    
    ###
    for i in range(k):
        DataPoints.pop(0)
    
    ###
    for iter in range(max_iter):
        
        indexingList = []
        for i in range(k):
            indexingList.append([])
            
        clustering(centroids_list ,DataPoints, indexingList )
        bol = update(centroids_list ,DataPoints, indexingList )
        
        if (bol == False):
            break
    return (centroids_list)


lst = kmeans(7 , 100 , r"C:\Users\weamm\Downloads\input_1.txt" , r"C:\Users\weamm\Downloads\out.txt")
print(lst)
"""
-4.3095,9.0182,5.3354
8.2052,-8.6995,-8.5803
9.7608,-5.7197,-7.2840
"""
