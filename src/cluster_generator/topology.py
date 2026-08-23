from itertools import combinations
import networkx as nx
from ase.data import atomic_numbers, covalent_radii


def get_metal_indices(atoms, metals):
    '''
    Return the ASE indices belonging to the user-defined cluster metals.
    '''
    return [atom.index for atom in atoms if atom.symbol in metals]


def get_metal_pair_distances(atoms, metal_indices):
    '''
    Calculate every unique metal-metal pair distance using MIC.
    '''
    pairs = []

    for i, j in combinations(metal_indices, 2):
        pairs.append({
            "i": i,
            "j": j,
            "symbol_i": atoms[i].symbol,
            "symbol_j": atoms[j].symbol,
            "distance": float(atoms.get_distance(i, j, mic=True)),
        })

    return pairs


def get_bond_cutoff(element_i, element_j, bond_cutoffs, bond_tolerance):
    '''
    Return the cutoff used to decide whether two metals are connected in the graph.
    '''
    pair_type = tuple(sorted((element_i, element_j)))

    if pair_type in bond_cutoffs:
        return bond_cutoffs[pair_type]

    radius_i = covalent_radii[atomic_numbers[element_i]]
    radius_j = covalent_radii[atomic_numbers[element_j]]
    return bond_tolerance * (radius_i + radius_j)


def build_metal_graph(atoms, metal_indices, pair_distances, bond_cutoffs, bond_tolerance):
    '''
    Build a NetworkX graph where metal atoms are nodes and connected metal pairs are edges.
    '''
    graph = nx.Graph()

    for index in metal_indices:
        graph.add_node(index, element=atoms[index].symbol)

    for pair in pair_distances:
        cutoff = get_bond_cutoff(pair["symbol_i"], pair["symbol_j"], bond_cutoffs, bond_tolerance)

        if pair["distance"] <= cutoff:
            graph.add_edge(pair["i"], pair["j"], distance=pair["distance"])

    return graph


def analyse_metal_graph(graph):
    '''
    Analyse connectivity, neighbours, coordination, edges and fully connected triangles.
    '''
    if graph.number_of_nodes() == 0:
        connected = False
    elif graph.number_of_nodes() == 1:
        connected = True
    else:
        connected = nx.is_connected(graph)

    neighbours = {node: list(graph.neighbors(node)) for node in graph.nodes}
    coordination = {node: graph.degree(node) for node in graph.nodes}
    edges = list(graph.edges)

    triangles = []
    for i, j, k in combinations(graph.nodes, 3):
        if graph.has_edge(i, j) and graph.has_edge(i, k) and graph.has_edge(j, k):
            triangles.append((i, j, k))

    return {
        "connected": connected,
        "neighbours": neighbours,
        "coordination": coordination,
        "edges": edges,
        "triangles": triangles,
    }


def get_growth_motifs(graph, graph_analysis):
    '''
    Return the metal centres, edges and triangles available for topological growth.
    '''
    return {
        "centres": list(graph.nodes),
        "edges": graph_analysis["edges"],
        "triangles": graph_analysis["triangles"],
    }


def analyse_structure(atoms, metals, bond_cutoffs, bond_tolerance):
    '''
    Run the current metal-framework analysis for one structure.
    '''
    metal_indices = get_metal_indices(atoms, metals)
    pair_distances = get_metal_pair_distances(atoms, metal_indices)
    metal_graph = build_metal_graph(atoms, metal_indices, pair_distances, bond_cutoffs, bond_tolerance)
    graph_analysis = analyse_metal_graph(metal_graph)
    growth_motifs = get_growth_motifs(metal_graph, graph_analysis)

    return {
        "atoms": atoms,
        "metal_indices": metal_indices,
        "n_metals": len(metal_indices),
        "pair_distances": pair_distances,
        "metal_graph": metal_graph,
        "connected": graph_analysis["connected"],
        "neighbours": graph_analysis["neighbours"],
        "coordination": graph_analysis["coordination"],
        "growth_centres": growth_motifs["centres"],
        "growth_edges": growth_motifs["edges"],
        "growth_triangles": growth_motifs["triangles"],
    }
