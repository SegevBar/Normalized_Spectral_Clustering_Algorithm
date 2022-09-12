import numpy as np
import spkmeansmodule
import sys


def main():
    '''
    Funcion: main()
    ----------------------------------------------------------------------------
    Params: input from CMD

    Action: Validate CMD input and runs Python program according goal

    Return: Output according to goal
    '''
    # validate only 4 CND args -> program, k, goal, filaname
    validateInput(len(sys.argv) != 4)

    # get k and validate it
    try:
        k = int(sys.argv[1])
    except:
        validateInput(False)
    validateInput('.' in str(k) or int(k) <= 0)

    # get goal and filename
    goal = sys.argv[2]
    filename = sys.argv[3]

    # run goal
    if goal == "spk":
        spk(k, filename)
    else:
        spkmeansmodule.goalsOtherThenSpk(goal, filename)

def spk(k, filename):
    '''
    Funcion: 
    ----------------------------------------------------------------------------
    Params: 

    Action: 

    Return: 
    '''
    # The function performs the normalizedSpectralClustering algorithm
    vectorsMatrix = spkmeansmodule.spkWithoutKmeans(k, filename)
    vectorsArray = np.array(vectorsMatrix)
    clusters = kmeanspp(vectorsArray)
    spkmeansmodule.kmeans(vectorsMatrix, clusters, np.size(vectorsArray, 1))


def kmeanspp(vectorsMatrix):
    '''
    Funcion: 
    ----------------------------------------------------------------------------
    Params: 

    Action: 

    Return: 
    '''
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

def printClusters(indices):
    '''
    Funcion: 
    ----------------------------------------------------------------------------
    Params: 

    Action: 

    Return: 
    '''
    # The function prints the indices
    for i in range(len(indices)):
        if i != len(indices) - 1:
            print(int(indices[i]), end=",")
        else:
            print(int(indices[i]))

def validateInput(bool):
    if bool:
        print("Invalid Input!")
        quit()

if __name__ == '__main__':
    main()

