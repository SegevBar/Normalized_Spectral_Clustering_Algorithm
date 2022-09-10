import numpy as np
import spkmeansmodule
import sys

'''
Funcion: 
-----------------------------------------------------------------------------
Params: 
Action: 
Return: 
'''
def main():
    # The function takes the input from argv and runs a function that
    # appropriate to the goal
    k = int(sys.argv[1])
    goal = sys.argv[2]
    filename = sys.argv[3]
    if goal == "spk":
        spk(k, filename)
    else:
        spkmeansmodule.goalsOtherThenSpk(goal, filename)

'''
Funcion: 
-----------------------------------------------------------------------------
Params: 
Action: 
Return: 
'''
def spk(k, filename):
    # The function performs the normalizedSpectralClustering algorithm
    vectorsMatrix = spkmeansmodule.spkWithoutKmeans(k, filename)
    VectorsArray = np.array(vectorsMatrix)
    clusters = kmeanspp(VectorsArray)
    spkmeansmodule.kmeans(vectorsMatrix, clusters, np.size(VectorsArray, 1))

'''
Funcion: 
-----------------------------------------------------------------------------
Params: 
Action: 
Return: 
'''
def kmeanspp(vectorsMatrix):
    # The function initializes K centroids for the K-means algorithm
    np.random.seed(0)
    len_column = np.size(vectorsMatrix, 0)
    len_row = np.size(vectorsMatrix, 1)

    #rearrange vectors for algorithm 
    vectorsMatrix = np.hstack(
        (np.arange(len_column).reshape(len_column, 1), vectorsMatrix))

    #prepare data structures
    clusters = np.zeros((len_row, len_row + 1))
    prob = [np.Infinity for i in range(len_column)]
    di = [np.Infinity for i in range(len_column)]
    lindex = np.random.choice(len_column)
    clusters[0] = vectorsMatrix[lindex, :]

    #main algorithm loop:
    for i in range(1, len_row):
        for j in range(len_column):
            temp = np.linalg.norm(vectorsMatrix[j, 1:] - clusters[i - 1][1:]) ** 2
            if temp < di[j]:
                di[j] = temp
        base = sum(di)
        for j in range(len_column):
            prob[j] = di[j] / base
        lindex = np.random.choice(len_column, p=prob)
        clusters[i] = vectorsMatrix[lindex, :]
    
    printClusters(clusters[:, 0])
    return clusters[:, 1:].tolist()

'''
Funcion: 
-----------------------------------------------------------------------------
Params: 
Action: 
Return: 
'''
def printClusters(indices):
    # The function prints the indices
    for i in range(len(indices)):
        if i != len(indices) - 1:
            print(int(indices[i]), end=",")
        else:
            print(int(indices[i]))

if __name__ == '__main__':
    main()

