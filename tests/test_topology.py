import networkx as nx
from cluster_generator.topology import analyse_metal_graph


def test_single_node_graph_is_connected():
    graph = nx.Graph()
    graph.add_node(0, element="Co")
    result = analyse_metal_graph(graph)

    assert result["connected"] is True
    assert result["coordination"][0] == 0


def test_triangle_is_detected():
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (0, 2)])
    result = analyse_metal_graph(graph)

    assert (0, 1, 2) in result["triangles"]
