def kmeans (k ,max_iter = 200, inputfile=None , outputfile=None):
    
    assert type(k) == int , "Invalid input."
    assert type(max_iter) == int , "Invalid input."
    
    file = open(outputfile)
    print(file.read())
    
    return 0
