#!/usr/bin/env python

"""
Author:
Student nr:
Script to:
"""

# Implement your functions here
global cycles, used_edges
cycles = []
used_edges = []

def init_globals():
    global cycles, used_edges
    cycles = []
    used_edges = []

def get_node_balances(graph):
    from_nodes = graph.keys()
    to_nodes   = graph.values()

    balances = {}

    for key in from_nodes:
        from_count = len(graph[key])

        to_count   = sum([1 for x in to_nodes if key in x])

        balances[key] = to_count - from_count

    return balances

def is_eulerian(balances):
    balances = [True if abs(x) < 1 else False for x in balances.values() ]
    if False in balances:   
        return False
    else: 
        return True

def find_next(to_nodes, from_node):
    for to_node in to_nodes:
        interaction = (from_node, to_node)
        if interaction not in used_edges:
            return interaction
    return None

def check_cycle(cycle):
    try:
        if cycle[0][0] == cycle[-1][1]:
            return True
        else:
            return False
    except:
        return False

def find_cycles(graph, from_node = None):
    split = []
    cycle = []
    if not from_node:
        from_node = graph.keys()[0]
    
    keep_reading = True
    while keep_reading:
        to_nodes = graph[from_node]
        #pick first approach
        if len(to_nodes) >= 1:
  
            interaction = find_next(to_nodes, from_node)
            if interaction:
                if len(to_nodes) > 1:
                    split.append(from_node)
                used_edges.append(interaction)
                cycle.append(interaction)
                
                from_node = interaction[1]
               
            else:
                keep_reading = False
        else:
            keep_reading = False


    if check_cycle(cycle):
        cycles.append(cycle)


    if len(split) > 0:
        find_cycles(graph, from_node = split[0])

def has_eulerian_path(balances):
    balances = [abs(x) for x in balances.values()]
    semi_c   = balances.count(1)
    if semi_c <= 2:
        return True    
    else:
        return False

def get_merge_idx(cycles):
    merge_idx = None
    try:
        merge_start = cycles[-1][0][0]
        

        for i in range(len(cycles[0])):
            if cycles[0][i][1] == merge_start:
                merge_idx = i
                break
    except IndexError:
        pass
    return merge_idx

def merge_cycles(cycles):
    cycles = sorted(cycles)
    merge_idx = get_merge_idx(cycles)
    if not merge_idx:
        return ""
    cycle_str = "".join(cycles[0][0])
    
    for i in range(1, len(cycles[0])):
        cycle_str += cycles[0][i][1]

        if i == merge_idx:
            for j in range(len(cycles[-1])):
                cycle_str += cycles[-1][j][1]

    return cycle_str

def get_eulerian_path(graph, balances):
    init_globals()

    start_key = None
    to_key = None
    for key,value in balances.items():
        if value   == 1:
            start_key = key
        elif value == -1:
            to_key    = key

    if to_key and start_key:
        graph[start_key].append(to_key)

    find_cycles(graph, start_key)
    raw_path = merge_cycles(cycles)
    return raw_path.split(start_key+to_key)[0]

if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = { 'A': ['B'],
                  'B': ['C'],
                  'I': ['H'],
                  'H': ['F'],
                  'F': ['G','E'],
                  'C': ['I','J'],
                  'G': ['A'],
                  'E': ['J'],
                  'J': ['F','D'],
                  'D': ['C']}


    graph_non = { 'A': ['B', 'C'],
                  'B': ['C'],
                  'I': ['H'],
                  'H': ['F'],
                  'F': ['G','E'],
                  'C': ['I','J'],
                  'G': ['A'],
                  'E': ['J'],
                  'J': ['F','D'],
                  'D': ['C']}

    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = { 5:  [6],
                     6:  [7],
                     10: [11],
                     11: [4],
                     4:  [5,3],
                     7:  [10,9],
                     3:  [9,1],
                     9:  [4,8],
                     8:  [7],
                     1:  [2],
                     2:  [3]}

    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']

    balances = get_node_balances(graph_non)



    if is_eulerian(balances):
        print("graph_822 is eulerian")
    else:
        print("graph_822 is not eulerian")

    if has_eulerian_path(balances):
        print("graph_822 has an eulerian path")
    else:
        print("graph_822 has no eulerian path")

    find_cycles(graph_non)
    print merge_cycles(cycles)

    print get_eulerian_path(graph_non, balances)