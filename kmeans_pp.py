def dist(point1, point2):
    sum = 0
    d = len(point1)
    for i in range(d):
        sum += (point1[i] - point2[i]) ** 2
    sum = sum ** 0.5
    return sum

def mindist(centroids_list, point, i):
    min = "inf"
    index = 0
    for j in range(i):
        distance = dist(centroids_list[j], point)
        if (min == "inf" or distance < min):
            min = distance
            index = j
    return index

