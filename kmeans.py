def dist(point1, point2):
    sum = 0
    d = point1.len()
    for i in range(d):
        sum += (point1[i] - point2[i]) ** 2
    sum = sum ** 0.5
    return sum


def mindist(centroids_list, point):
    min = "inf"
    index = 0
    for j in range(centroids_list.len()):
        distance = dist(centroids_list[j], point)
        if (distance < min or min == "inf"):
            min = distance
            index = j
    return index


def clustering(centroids_list, DataPoints, indexingList):
    for i in range(DataPoints.len()):
        x = mindist(centroids_list, DataPoints[i])
        indexingList[x].append(i)


def update(centroids_list, DataPoints, indexingList):
    boolean = True
    for i in range(centroids_list.len()):
        old_avg = centroids_list[i]
        centroids_list[i] = Avg(indexingList[i], DataPoints)
        distance = dist(centroids_list, old_avg)
        if (distance > 0.001):
            boolean = False
    return boolean


def Avg(lst, DataPoints):
    avg_point = []
    d = DataPoints[0].len()
    for i in range(d):
        sum = 0
        for item in lst:
            sum += DataPoints[item][i]
        sum = sum / lst.len()
        avg_point.append(sum)
    return avg_point


def kmeans(k, max_iter=200, inputfile=None, outputfile=None):
    assert type(k) == int, "Invalid input."
    assert type(max_iter) == int, "Invalid input."

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

    for iter in range(max_iter):

        indexingList = []
        for i in range(k):
            indexingList.append([])

        clustering(centroids_list, DataPoints, indexingList)
        bol = update(centroids_list, DataPoints, indexingList)
        if (bol == False):
            break

    return (centroids_list)


kmeans(2, 100, r"C:\Users\weamm\Downloads\input_1.txt", r"C:\Users\weamm\Downloads\out.txt")

"""
-4.3095,9.0182,5.3354
8.2052,-8.6995,-8.5803
9.7608,-5.7197,-7.2840
"""

