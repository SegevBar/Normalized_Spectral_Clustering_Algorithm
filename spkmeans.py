import numpy as np
import sys
import spkmeansmodule


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
    validateInput('.' in str(k) or int(k) < 0)

    # get goal and filename
    goal = sys.argv[2]
    filename = sys.argv[3]

    # run goal
    if goal == "spk":
        spk(k, filename)
    else:
        spkmeansmodule.runGoalOfCProgram(goal, filename)


def spk(k, filename):
    '''
    Funcion: spk(k, filename)
    ----------------------------------------------------------------------------
    Params: k, input file

    Action: Performs the Normalized Spectral Clustering algorithm

    Return: Prints indexes of initial centroids and calculated centroids
    '''
    # The function performs the normalizedSpectralClustering algorithm
    vectorsMatrix = spkmeansmodule.getPythonNormalizedKEigenvectorsMatrix\
        (k, filename)
    
    # init parameters and run kmeans++
    vectorsArray = np.array(vectorsMatrix)
    n = np.size(vectorsArray, 0)
    vectorDim = np.size(vectorsArray, 1)

    centroids = np.array([[0.0 for i in range(vectorDim)] \
        for i in range(vectorDim)])
    chosenKCentroidsIndexs = kmeanspp(vectorsArray, centroids, n, vectorDim)
    
    # init parameters and run kmeans from C program
    kmeansArgs = getArgsForKmeans(vectorsArray, centroids)
    kmeansCentroids = np.array(spkmeansmodule.runKmeansFromCProgram(\
        kmeansArgs[0],kmeansArgs[1], n, vectorDim))
    print(kmeansCentroids)
    printClusters(chosenKCentroidsIndexs, kmeansCentroids)



def kmeanspp(vectorsMatrix, centroids, n, k):
    '''
    Funcion: kmeanspp(vectorsMatrix, centroids)
    ----------------------------------------------------------------------------
    Params: numpy matrix of c program normalized k eigenvectors matrix (T)

    Action: Init paramameters and run kmeans++ algorithm
            Updates centroids list of lists

    Return: List of chosen centroids indexes
    '''
    np.random.seed(0)

    # prepare data structures
    distances = [-1.0 for i in range(n)]
    probs = [0.0 for i in range(n)]
    centroidsLoc = [0 for i in range(k)]

    # choose a random datapoint to be the first centroid
    rendIndex = np.random.choice(n)
    centroids[0] = np.ndarray.copy(vectorsMatrix[rendIndex])
    centroidsLoc[0] = rendIndex

    # algorithm loop
    for i in range(1, k):
        minDistance(vectorsMatrix, centroids, distances, i)
        calcProbability(vectorsMatrix, distances, probs)
        rendIndex = np.random.choice(n, p = probs)
        centroidsLoc[i] = rendIndex
        centroids[i] = np.ndarray.copy(vectorsMatrix[rendIndex])

    return centroidsLoc


def minDistance(vectorsMatrix, centroids, distances, D):
    '''
    Funcion: minDistance(vectorsMatrix, centroids, distances, D)
    ----------------------------------------------------------------------------
    Params: data points, centroids vectors list, calculated distances and 
            current centroid index

    Action: Calculates Dl = min(xl-uj)^2 
            Updates distances list accordingly

    Return: None
    '''
    for i in range(len(vectorsMatrix)):
        curDistance = pow(np.linalg.norm(vectorsMatrix[i] - centroids[D-1]), 2)
        
        if curDistance < distances[i] or distances[i] == -1.0:
            distances[i] = curDistance


def calcProbability(vectorsMatrix, distances, probs):
    '''
    Funcion: calcProbability(vectorsMatrix, distances, probs)
    ----------------------------------------------------------------------------
    Params: data points, calculated distances, calculated probabilities
            
    Action: Calculates P(xl) = Dl/sum(Dm) m=i,..,N
            Updates probs list accordingly

    Return: None
    '''
    sum = 0.0
    for dist in distances:
        sum += dist
    for i in range(len(vectorsMatrix)):
        probs[i] = distances[i] / sum


def getArgsForKmeans(vectorsMatrix, centroids):
    '''
    Funcion: getArgsForKmeans(vectorsMatrix, centroids)
    ----------------------------------------------------------------------------
    Params: data points, centroids vectors list

    Action: Converts vectors matrix and centroids matrix to 1D list

    Return: tuple of args (dataPoints1D, centroids1D, n, dim)
    '''
    #set arguments for c extension use
    dim = len(centroids[0])
    n = len(vectorsMatrix)

    #transform to 1D list of data points for the c extension
    dataPoints1D = []
    for vec in vectorsMatrix:
        for i in range(dim):
            dataPoints1D.append(vec[i])

    #transform to 1D list of centroids for the c extension
    centroids1D = []
    for centroid in centroids:
        for i in range(dim):
            centroids1D.append(centroid[i])
    
    return (dataPoints1D, centroids1D, n, dim)
    

def printClusters(chosenKCentroidsIndexs, kmeansCentroids):
    '''
    Funcion: printClusters(indices)
    ----------------------------------------------------------------------------
    Params: List of kmeans++ centroid indexs, centroids vectors

    Action: Prints centroid indexs and kmeans centroids

    Return: None
    '''
    print(chosenKCentroidsIndexs, sep=",")
    for centroid in kmeansCentroids:
        print(*["{:.4f}".format(num) for num in centroid], sep=",")


def validateInput(bool):
    '''
    Funcion: validateInput(bool)
    ----------------------------------------------------------------------------
    Params: boolean

    Action: Abort program if the input isn't valid

    Return: None
    '''
    if bool:
        print("Invalid Input!")
        quit()


if __name__ == '__main__':
    main()

