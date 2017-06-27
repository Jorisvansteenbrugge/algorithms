#!/usr/bin/env python

"""
Author: Joris van Steenbrugge
Student number: 950416
Implementation of the k-means clustering algorithm
"""

from __future__ import division
from math import pow, sqrt
from operator import mul
import random   


def get_dim_mean(points, idx):
    return sum([point[idx] for point in points]) / len(points)

def get_centroid(points, cluster):
    points = [point for point in points if point[-1] == cluster]
    return [get_dim_mean(points, idx) for idx in range(len(points[0]) - 1) ]



def euclidean_distance(x, y):
    power = 2
    try:
        return sqrt(sum([pow((x[i] - y[i]), power) for i in range(len(x))]))
    except TypeError:
        return None
    

def initialize_cluster(point, k_range):
    point[-1] = random.choice(k_range)
    return point
    
def csv_parser(lines):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """ 

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(map(float, items[1:])) #first item is the label
        except ValueError: #must be the header
            continue
    return data_points

def structure_points(points):
    return  [x + [None] for x in points]

def re_assign_point(point, centroids):
    """Return datapoints with re-assigned clusters.

        Keyword Arguments:
            point     -- list, datapoint containing the measurement for each 
                               dimension according to [x, y, ..., cluster]. 
            centroids -- list, containing the central points of each cluster
                               (based on averages). 

        The euclidean distance between the datapoint and all centroids is 
        calculated. The datapoint is then (re-)assigned to the closest centroid's 
        cluster.
    """
    point_dims = len(point) - 1 # last one is cluster number
    distances  = [euclidean_distance(point[0: point_dims], centroid) for \
                                centroid in centroids]
    closest    =  distances.index(min(distances))

   

    point[-1]  = closest
    return point

def random_centroid(points):
    """Return a random point in the dataset to serve as centroid.

        Keyword Arguments:
            points -- list, each entry is a list with datapoint 
                            measurements. Each datapoint is according
                            to: [x, y, ..., cluster] where the number of 
                            measurements is variable. The last value is 
                            always the clustering number. 
    """
    centroid = random.choice(points)
    return centroid[0: len(centroid) - 1]


def centroid_change(c, prev_c, cutoff):
    """Returns True if the difference is lower than the cutoff.
        
        Keyword Arguments:
            c      -- list, containing the central point of a cluster
                            (based on averages). 
            prev_c -- list, containing the central points of a cluster
                            (based on averages) for the centroid in the
                            previous iteration.
            cutoff -- float, A cutoff value for approximate convergence. 
                             If the percentage change of the centroid 
                             between iterations is lower than the cutoff 
                             the algorithm is artificially forced to 
                             "converge".

        The difference between c and prev_c is calculated as the percentage
        difference for each dimension. If the value in all dimensions is lower
        than the cutoff, True is returned.
    """
    abs_dim_change = [abs(c[i] - prev_c[i])         for i in range(len(c))]
    percent_change = [abs_dim_change[i] / prev_c[i] for i in range(len(c))]

    check = [True if percent < cutoff else False 
                  for percent in percent_change]

    if False not in check: # Stop when this happens
        return True
    else:
        return False


def calc_centroid_change(centroids, prev_centroids, cutoff):
    """Return Trues if centroid change is lower than the cutoff else Falses
        
        Keyword Arguments:
            centroids      -- list, containing the central points of each cluster
                                    (based on averages). 
            prev_centroids -- list, containing the central points of each cluster
                                    (based on averages) for the centroids in the
                                    previous iteration.
            cutoff         -- float, A cutoff value for approximate convergence. 
                                     If the percentage change of the centroid 
                                     between iterations is lower than the cutoff 
                                     the algorithm is artificially forced to 
                                     "converge".

        Each value in the return is a list with a (combination) of Trues and 
        Falses. The length of the nested lists is equal to the dimensions
        of the datapoints. If all values in an individual nested list are True
        it means that the difference with the previous centroid was lower than
        the cutoff.

    """
    return [centroid_change(centroids[i], prev_centroids[i],
                        cutoff) for i in range(len(centroids))]

def Kmeans(points, k = 2, cutoff = 0.02):
    """Returns clustered data-points using the Kmeans clustering algorithm.

        Keyword Arguments:
            points -- list, each entry is a list with datapoint 
                                    measurements. Each datapoint is according
                                    to: [x, y, ..., cluster] where the number of 
                                    measurements is variable. The last value is 
                                    always the clustering number.
            k      -- int, The number of clusters to be created.
            cutoff -- float, a cutoff value for approximate convergence. If the
                        percentage change of the centroid between iterations is
                        lower than the cutoff the algorithm is artificially 
                        forced to "converge".
    """
    clusters = range(k)

    # Each cluster is assigned a centroid at the start, this does not prevent,
    # the cluster from getting emptied later.
    centroids = [random_centroid(points)        for cluster in clusters]


    iter_times      = 0
    prev_centroids  = centroids
    while True:
        points      = [re_assign_point(point, centroids) for point   in points]

        try:
            centroids  = [get_centroid(points, cluster)  for cluster in clusters]
        # This happens because at first, none of the data-points are assigned to
        # a cluster.
        except IndexError:
            pass
       
        iter_times += 1


        # Stop the while loop if the algorithm has converged
        if centroids == prev_centroids:    
            print("K-means converged in {} iterations".format(iter_times))
            return points

        # Let the algorithm approximately converge if the precentage change
        # is lower than a cutoff.
        else:
            check = calc_centroid_change(centroids, prev_centroids, cutoff)
            if False not in check:
                print("K-means approximately converged in {} iterations".format(
                            iter_times))
                return points

        
        prev_centroids = centroids


def prepare_datapoints(file_name):
    """Return structured data-points loaded from a file

        Keyword Arguments:
            file_name -- str, file name, containing a datapoint on each line
                            with an x number of dimensions.

        This function serves as a simple wrapper.
    """
    points    = csv_parser(open(file_name))
    points    = structure_points(points)
    return points


def pretty_print(clustered_points):
    """Print data points in a human readable way.

        Keyword Arguments:
            clustered_points -- list, each entry is a list with datapoint 
                                    measurements. Each datapoint is according
                                    to: [x, y, ..., cluster] where the number of 
                                    measurements is variable. The last value is 
                                    always the clustering number.

    """
    clustered_points = ["\t".join(map(str, x)) for x in clustered_points]
    print "\n".join(clustered_points)

if __name__  == "__main__":
    points = prepare_datapoints("/home/joris/Downloads/datasets/2dtest.csv")
    clustered_points = Kmeans(points, 5, cutoff = 0.02)
    #pretty_print(clustered_points)

    for i in [2,3,4,5,6]:
        points = prepare_datapoints("/home/joris/Downloads/datasets/LargeSet_1.csv")
        clustered_points = Kmeans(points, i, cutoff = 0.02)
     #   pretty_print(clustered_points)

    points = prepare_datapoints("/home/joris/Downloads/datasets/LargeSet_2.csv")
    clustered_points = Kmeans(points, 2, cutoff = 0.02)
    #pretty_print(clustered_points)
 