#!/usr/bin/env python
"""
Author: Joris van Steenbrugge
Student number: 950416
Implementation of the k-means clustering algorithm

Hints:
- write a function to obtain Euclidean distance between two points.
- write a function to initialize centroids by randomly selecting points 
    from the initial set of points. You can use the random.sample() method
- write a function to find the closest centroid to each point and assign 
    the points to the clusters.     
- write a function to calculate the centroids given the clusters
"""
from __future__ import division
from math import pow, sqrt
import random   


def get_dim_mean(points, idx):
    return sum([point[idx] for point in points]) / len(points)

def get_centroid(points, cluster):
    points = [point for point in points if point[-1] == cluster]
    
    return [get_dim_mean(points, idx) for idx in range(len(points[0]) - 1) ]


def euclidean_distance(x, y):
    power = 2
    return sqrt(sum([pow((x[i] - y[i]), power) for i in range(len(x))]))
    

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
    return [x + [None] for x in points]

def re_assign_point(point, centroids):
    point_dims = len(point) - 1 # last one is cluster number
    distances  = [euclidean_distance(point[0:point_dims], centroid) for \
                                centroid in centroids]
    closest    =  distances.index(max(distances))
    point[-1] = closest
    return point


def Kmeans(points, k = 2):
    clusters = range(k)
    points = [initialize_cluster(point, clusters) for point in points]

    c = 0

    prev_points = []
    while True:
        centroids = [get_centroid(points, cluster) for cluster in clusters]
        points = [re_assign_point(point, centroids) for point in points]

        if points == prev_points:
            return points

        prev_points = points

        print c
        c+= 1


if __name__ == "__main__":
    test_file = "/home/joris/Downloads/datasets/2dtest.csv"
    points = csv_parser(open(test_file))
    points = structure_points(points)

    clustered_points = Kmeans(points)
    print clustered_points