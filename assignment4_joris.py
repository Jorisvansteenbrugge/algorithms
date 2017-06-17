#!/usr/bin/env python

"""
Author:
Student nr:
Script to:
"""

# Implement your functions here

def find_next(to_nodes, from_node,used_edges):
    for to_node in to_nodes:
        interaction = (from_node, to_node)
        if interaction not in used_edges:
            from_node = to_node
            return interaction
    return None


def is_eulerian(graph):
    used_edges = []
    from_node = graph.keys()[0]

    keep_reading = True
    while keep_reading:
        to_nodes = graph[from_node]
        #pick first approach
        if len(to_nodes) >= 1:
            interaction = find_next(to_nodes, from_node, used_edges)
            if interaction:
                used_edges.append(interaction)
                from_node = interaction[1]
            else:
                keep_reading = False
        else:
            keep_reading = False
    #check if first and last nodes are equal     
    if used_edges[0][0] == used_edges[-1][1]:
       return True
    else:
        return False

if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','J'],'G':['A'],'E':['J'],'J':['F','D'],'D':['C']}
    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = {5:[6],6:[7],10:[11],11:[4],4:[5,3],\
        7:[10,9],3:[9,1],9:[4,8],8:[7],1:[2], 2:[3]}
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']



    if is_eulerian(graph_822):
        print "graph_822 is eulerian"
    else:
        print "graph_822 is not eulerian"

    if is_eulerian(bigger_graph):
        print "bigger_graph is eulerian"
    else:
        print "bigger_graph is not eulerian"
    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file) 
