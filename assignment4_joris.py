#!/usr/bin/env python

"""
Author: Joris van Steenbrugge
Student nr: 950416798110
Script to: Find a Eulerian path or cycle in a graph
"""

global cycles, used_edges, split
cycles = []
used_edges = []
split = []

def init_globals():
    """Re-initializes global variables to empty lists
    """
    global cycles, used_edges, split
    cycles = []
    used_edges = []
    split = []

def get_node_balances(graph):
    """Returns the balances for each node in a graph.

        Keyword Arguments:
            graph -- dictionary, {node_label: [connectionA, connectionB]}.
    """
    from_nodes = graph.keys()
    to_nodes   = graph.values()

    balances = {}

    for key in from_nodes:
        from_count = len(graph[key])

        to_count   = sum([1 for x in to_nodes if key in x])

        balances[key] = to_count - from_count

    return balances

def is_eulerian(balances):
    """Returns if the graph is Eulerian based on the balances.

        Keyword Arguments: 
            balances -- dictionary, {node:balance}

        A graph is considered Eulerian if all nodes are balanced.

    """
    balances = [True if abs(x) < 1 else False for x in balances.values() ]
    if False in balances:   
        return False
    else: 
        return True

def find_next(to_nodes, from_node):
    """Return the next node

        Keyword Arguments:
            to_nodes  -- string, 
            from_node -- 
    """
    for to_node in to_nodes:
        interaction = (from_node, to_node)
        if interaction not in used_edges:
            return interaction
    return None

def check_cycle(cycle):
    """Returns a boolean corresponding to if the path is a cycle or not.

        Keyword Arguments:
            cycles -- list, containing possible cycles found in the graph
                        e.g. cycles: [cycle1, cycle2]
                            cycle1: [('A', 'T'), ('T', 'C')] 
    """
    try:
        if cycle[0][0] == cycle[-1][1]:
            return True
        else:
            return False
    except:
        return False

def find_cycles(graph, from_node = None):
    """Find cycles in a graph containing a Eulerian cycle and/or path.

        Keyword Arguments:
            graph     -- dictionary, {node_label: [connectionA, connectionB]}.
            from_node -- string, the node to start the cycle/path with.

        Returns:
            split -- A variable containing the nodes that connect with multiple
                     other nodes.
    """
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

    return split

def has_eulerian_path(balances):
    """Returns if the graph has a Eulerian path based on the balances.

        Keyword Arguments: 
            balances -- dictionary, {node:balance}

        A graph has a Eulerian path if at most two nodes are semi-balanced
        and all other nodes are balanced.
            
    """
    balances = [abs(x) for x in balances.values()]
    semi_c   = balances.count(1)
    if semi_c <= 2:
        return True    
    else:
        return False

def get_merge_idx(cycles):
    """Returns the index of an overlapping point in the graph to merge 2 cycles
        
        Keyword Arguments:
            cycles -- list, containing possible cycles found in the graph
                        e.g. cycles: [cycle1, cycle2]
                            cycle1: [('A', 'T'), ('T', 'C')] 
    """
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
    """Returns a merged Eucledian- cycle or path.

        Keyword Arguments:
            cycles -- list, containing possible cycles found in the graph
                        e.g. cycles: [cycle1, cycle2]
                            cycle1: [('A', 'T'), ('T', 'C')] 

    """
    cycles = sorted(cycles)

    if len(cycles) == 1:
        cycle_out = list(cycles[0][0][0])
        cycle_out += [cycle[1] for cycle in cycles[0]]
        return cycle_out

    merge_idx = get_merge_idx(cycles)
    
    if not merge_idx:
        return ""

    cycle_out = list(cycles[0][0])
    

    for i in range(1, len(cycles[0])):
        cycle_out.append(cycles[0][i][1])

        if i == merge_idx:
            for j in range(len(cycles[1])):
                cycle_out.append(cycles[1][j][1])
                

    return cycle_out

def get_eulerian_path(graph, balances):
    """Returns a Eulerian path in a graph.

        Keyword Arguments:
            graph -- dictionary, {node_label: [connectionA, connectionB]}

        The graph does not have to be Eulerian, at max 2 nodes 
        can be semi-balanced.
    """
    init_globals()

    start_key = None
    to_key = None

    # Find the semi-balanced node keys
    for key,value in balances.items():
        if value   == 1:
            start_key = key
        elif value == -1:
            to_key    = key

    # If there are semi-balanced node keys, an edge is artificially added
    # between these two nodes.
    was_eulerian = True
    if to_key and start_key:
        was_eulerian = False
        graph[start_key].append(to_key)

    splits = find_cycles(graph, start_key)

    # This is probably not necessary.
    for node in splits:
        find_cycles(graph, node)


    raw_path = merge_cycles(cycles)
    
    if was_eulerian: # Then the path will be equal to the cycle
        return "".join(raw_path)
    else: # We have to remove the artificially added node 
        path = "".join(raw_path)
        return "".join(path.split(start_key+to_key)[::-1])

def check_graph(graph):
    """Returns a graph with integers as key/values as a string graph
        
        Keyword Arguments:
            graph -- dictionary, {node_label: [connectionA, connectionB]}
    """
    string_graph = {}
    for key in graph.keys():
        string_graph[str(key)] = map(str, graph[key])
    return string_graph

def process_graph(graph):
    """Wrapper function to find a Eulcidean cycle or path in a graph

        Keyword Arguments:
            graph -- dictionary, {node_label: [connectionA, connectionB]}

        It was not until later that I found out that my script did not
        handle integer graphs well out of the box, this function is some kind
        of quick and dirty solution.
    """
    init_globals()

    # If the graph contains integers this is triggered.
    # See the check_graph() function description for more information
    if type(graph.keys()[0]) == int:
        graph = check_graph(graph)

    
    balances = get_node_balances(graph)


    if is_eulerian(balances):
        splits = find_cycles(graph, sorted(graph.keys())[0])
        for node in splits:
            find_cycles(graph, node)

        print("graph is eulerian:\n{}\n".format(
            "->".join(merge_cycles(cycles))))


    else:
        print("graph is not eulerian")

    if has_eulerian_path(balances):
        print("graph has an eulerian path:\n{}\n".format(
            "->".join(get_eulerian_path(graph, balances))))
    else:
        print("graph has no eulerian path")

def graph_from_spectrum(spectrum):
    """Returns a graph based on a spectrum with kmers.

        Keyword Arguments:
            spectrum -- list of strings, containing kmers.

        Each node will consist of the fist or second half of the
        kmer (the length of each node sequence is l-1).
    """
    graph = {}

    # init empty
    for kmer in spectrum:
        l = len(kmer)
        k1 = kmer[0: l - 1]
        k2 = kmer[l - 2:  ]

        if k1 not in graph:
            graph[ k1 ] = [k2]
        else:
            graph[ k1 ].append(k2)
        if k2 not in graph:
            graph[ k2 ] = []
        else:
            pass
        


    return graph
    
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

    process_graph(graph_822)
    
    graph_spectrum = graph_from_spectrum(s)
 
    print("The graph based on the spectrum:")
    for k,v in graph_spectrum.items():
        print(k,v)

    process_graph(graph_spectrum)


    process_graph(bigger_graph)