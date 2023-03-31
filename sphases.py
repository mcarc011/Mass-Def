#%%
import matplotlib
import itertools
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os


def tako_sample_arrays() -> tuple:
    """Sample arrays.

    Returns:
        tuple: Sample arrays.
    """
    black_line_mtx = np.array([[0,2,0,0],[0,0,0,0],[0,0,0,0],[2,0,0,0]])
    red_line_mtx = np.array([[0,0,2,0],[0,0,0,0],[2,0,0,0],[0,0,0,0]])
    b2 = np.array([[0,2,3,4],[5,6,7,8],[0,0,0,0],[2,0,0,0]])
    return (black_line_mtx, red_line_mtx, b2)

def mkdir_if_dne(path: str):
    """Make directory if it does not exist.

    Args:
        path (str): Path to directory.
    """
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    return

def get_edges(input_mtx: np.array) -> dict:
    """Determine what edges are connect, based on matrix elements

    Args:
        input_mtx (np.array): Matrix describing connected elements.

    Returns:
        dict: Ordered dictionary with keys of nodes in connected order and values of "weights"/"number of arrows".
    """
    edge_labels_dict = {}
    for i in range(len(input_mtx)):
        for j in range(len(input_mtx[i])):
            if i != j:
                number_of_arrows = input_mtx[i][j]
                if number_of_arrows != 0:
                    start_node = i+1
                    end_node = j+1
                    edge_labels_dict[(start_node, end_node)] = int(number_of_arrows)
    return edge_labels_dict

def set_plotting_style(style="default"):
    """(Unused function.) Set plot style.  Does not work with nx.

    Args:
        style (str, optional): Style for plots. Defaults to "default".
    """
    plt.style.use(style)
    return

def draw_edges_and_labels(G: nx.DiGraph, pos: dict, edge_labels_dict: dict, color: str, font_size=12, node_color="yellow", node_size=500, label_pos=0.1, **kwargs):
    """Draw the connected nodes.

    Args:
        G (nx.DiGraph): networkx graph.
        pos (dict): Positions of nodes.
        edge_labels_dict (dict): Ordered dictionary with keys of nodes in connected order and values of "weights"/"number of arrows".
        color (str): Edge and label's font color.
        font_size (int, optional): Font size of label. Defaults to 12.
        node_color (str, optional): Color of nodes. Defaults to "yellow".
        node_size (int, optional): Node size. Defaults to 500.
        label_pos (float, optional): Position of edge label along edge (0=head, 0.5=center, 1=tail). Defaults to 0.1.
    """
    nx.draw(
        G, pos, edge_color=color, width=1, linewidths=0.5,
        node_size=node_size, node_color=node_color, alpha=0.9,
        labels={node: node for node in G.nodes()},
        **kwargs
    )
    if (color=="red") | (color=="r"):
        nx.draw_networkx_edges(G, pos, arrowsize=1, edge_color=color)
        nx.draw_networkx_edge_labels(
            G, pos, label_pos=0.5,
            edge_labels=edge_labels_dict,
            font_color=color,
            font_size=font_size,
            **kwargs
        )
    else:
        nx.draw_networkx_edges(G, pos, arrowsize=20, edge_color=color)
        nx.draw_networkx_edge_labels(
            G, pos, label_pos=label_pos,
            edge_labels=edge_labels_dict,
            font_color=color,
            font_size=font_size,
            **kwargs
        )
    return

def add_edges_to_graph(G: nx.DiGraph, edge_labels_dict: dict, color: str):
    """Add edges to connected nodes on graph.

    Args:
        G (nx.DiGraph): networkx graph.
        edge_labels_dict (dict): Ordered dictionary with keys of nodes in connected order and values of "weights"/"number of arrows".
        color (str): Color of edges.
    """
    list_for_add_edges_from = list(edge_labels_dict.keys())
    G.add_edges_from(list_for_add_edges_from, color=color)
    return

def generate_quiver_diagram(input_mtx: np.array, black_edge_labels_dict: dict, black_color: str,  figsize=(5,5), fname='test'):
    """Generate quiver diagrams.  Connected nodes are determined by input matricies. 

    Args:
        input_mtx (np.array): Matrix of nodes connected with black lines.
        black_edge_labels_dict (dict): Black edge ordered dictionary with keys of nodes in connected order and values of "weights"/"number of arrows".
        black_color (str): Color used for plot black edges and labels.
        red_edge_labels_dict (dict): Red edge ordered dictionary with keys of nodes in connected order and values of "weights"/"number of arrows".
        red_color (str): Color used for plot black edges and labels.
        figsize (tuple): Size of output graphic. Defaults to (12,8).
        savefig (str or False, optional): Path for saved figure.  When set to False, figure is shown but not saved. Defaults to False.
    """
    fig = plt.figure(figsize=figsize)

    # set_plotting_style("ggplot") # doesnt work with nx
    Gk = nx.DiGraph() # for black nodes / edges
    Gr = nx.DiGraph() # for red nodes / edges
    
    # NOTE: shows connected nodes...ALSO SHOWS UNCONNECTED NODES (if they exist)
    nodes = np.arange(1, len(input_mtx)+1).tolist()
    Gk.add_nodes_from(nodes)

    # # NOTE: ONLY SHOWS CONNECTED NODES
    # list_for_add_edges_from = list(edge_labels_dict.keys())
    # G.add_edges_from(list_for_add_edges_from)

    # add black edges
    add_edges_to_graph(Gk, black_edge_labels_dict, black_color)

    # choose node locations / layout
    pos = nx.circular_layout(Gk)

    # add edges and labels
    draw_edges_and_labels(Gk, pos, black_edge_labels_dict, black_color, label_pos=0.3)

    plt.savefig(fname)
    plt.show()
    plt.close()
    return


def plot_web(X: np.array, fn: str):
    black_line_mtx =X    
    k_edge_labels_dict = get_edges(black_line_mtx)
    generate_quiver_diagram(black_line_mtx, k_edge_labels_dict, "k", fname=fn)
    return

text = ''
for pi,phase in enumerate(dweb):
    plot_web(phase[0],fn=str(pi))
    text+= '''
|<img src="./''' 
    text += str(pi)
    text+='''.png" width="200" height="200"> |
|---|
|Phase '''
    text+= str(pi+1)
    text+= '|'
    Wstr = '+'.join(phase[-1])
    Wstr = Wstr.replace(',',' ')
    text+= '\n\nSuperpotential: ' + str(pi+1) +'\n' + Wstr+'\n\n'

f=open('README2.md','w')
f.write(text)
f.close()
# %%
