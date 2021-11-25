def dist(point1, point2):
    sum = 0
    d = len(point1)
    for i in range(d):
        sum += (point1[i] - point2[i]) ** 2
    sum = sum ** 0.5
    return sum


def mindist(centroids_list, point):
    min = "inf"
    index = 0
    for j in range(len(centroids_list)):
        distance = dist(centroids_list[j], point)
        if (min == "inf" or distance < min):
            min = distance
            index = j
    return index


def clustering(centroids_list, DataPoints, indexingList):
    for i in range(len(DataPoints)):
        x = mindist(centroids_list, DataPoints[i])
        indexingList[x].append(i)


def update(centroids_list, DataPoints, indexingList):
    boolean = True
    for i in range(len(centroids_list)):
        old_avg = centroids_list[i]
        centroids_list[i] = Avg(indexingList[i], DataPoints)
        distance = dist(centroids_list[i], old_avg)
        if (distance > 0.001):
            boolean = False

    return boolean


def Avg(lst, DataPoints):
    avg_point = []
    d = len(DataPoints[0])
    for i in range(d):
        sum = 0
        for item in lst:
            sum += DataPoints[item][i]
        sum = sum / len(lst)
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

        if (bol == True):
            break

    for i in range(len(centroids_list)):
        for j in range(len(centroids_list[0])):
            a = centroids_list[i][j]
            centroids_list[i][j] = "%.4f" %round(a, 4)

    return (centroids_list)


lst = kmeans(15, 100, r"C:\Users\nadsa\Documents\input_3.txt")
f = open(r"C:\Users\nadsa\Documents\result_1.txt", "w")
for l in lst:
    for i in range(len(l)):
        f.write(str(l[i]))
        if(i != len(l)-1):
            f.write(",")
    f.write("\n")


"""
-4.3095,9.0182,5.3354
8.2052,-8.6995,-8.5803
9.7608,-5.7197,-7.2840
"""
