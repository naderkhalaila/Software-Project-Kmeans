def kmeans (k ,max_iter = 200, inputfile=None , outputfile=None):
    
    assert type(k) == int , "Invalid input."
    assert type(max_iter) == int , "Invalid input."
    
    file = open(inputfile)
    content = file.read()
    lst = content.splitlines()
    for i in range(len(lst)):
        lst[i] = lst[i].split(",")
        for j in range(len(lst[i])):
            lst[i][j] = float(lst[i][j])
    
    return 0



