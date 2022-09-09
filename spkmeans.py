import numpy as np
import spkmeansmodule
import sys

def kmeanspp(dataPoints):
    # The function initializes K centroids for the K-means algorithm
    np.random.seed(0)
    len_column = np.size(dataPoints, 0)
    len_row = np.size(dataPoints, 1)

    #rearrange dapa points for algorithm 
    dataPoints = np.hstack(
        (np.arange(len_column).reshape(len_column, 1), dataPoints))

    #prepare data structures
    clusters = np.zeros((len_row, len_row + 1))
    prob = [np.Infinity for i in range(len_column)]
    di = [np.Infinity for i in range(len_column)]
    lindex = np.random.choice(len_column)
    clusters[0] = dataPoints[lindex, :]

    #main algorithm loop:
    for i in range(1, len_row):
        for j in range(len_column):
            temp = np.linalg.norm(dataPoints[j, 1:] - clusters[i - 1][1:]) ** 2
            if temp < di[j]:
                di[j] = temp
        base = sum(di)
        for j in range(len_column):
            prob[j] = di[j] / base
        lindex = np.random.choice(len_column, p=prob)
        clusters[i] = dataPoints[lindex, :]
    
    printClusters(clusters[:, 0])
    return clusters[:, 1:].tolist()


def printClusters(indices):
    # The function prints the indices
    for i in range(len(indices)):
        if i != len(indices) - 1:
            print(int(indices[i]), end=",")
        else:
            print(int(indices[i]))


def main():
    # The function takes the input from argv and runs a function that
    # appropriate to the goal
    k = int(sys.argv[1])
    goal = sys.argv[2]
    filename = sys.argv[3]
    if goal == "spk":
        spk(k, filename)
    else:
        spkmeansmodule.goalsOtherThenSpk(k, goal, filename)


def spk(k, filename):
    # The function performs the normalizedSpectralClustering algorithm
    dataPoints = spkmeansmodule.spkWithoutKmeans(k, filename)
    dataPointsArray = np.array(dataPoints)
    clusters = kmeanspp(dataPointsArray)
    spkmeansmodule.kmeans(dataPoints, clusters, np.size(dataPointsArray, 1))


if __name__ == '__main__':
    main()

