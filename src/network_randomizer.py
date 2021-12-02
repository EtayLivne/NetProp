from random import choice
from multiprocessing import Queue, Pool
import os.path as os_path
import common.network_loaders as network_loaders
import networkx as nx

def generate_rank_equivalent_random_network(network, edge_switch_factor):
    switches = edge_switch_factor * len(list(network.edges))
    switch_types = ["place", "target", "source"]
    edge_source_sample = []
    for node in network.nodes:
        edge_source_sample.extend([node] * network.degree(node))

    while switches:
        switch_type = choice(switch_types)
        source_1 = choice(edge_source_sample)
        target_1 = choice(list(network.neighbors(source_1)))
        edge_1 = (source_1, target_1)
        edge_1_attrs = dict(network.edges[edge_1])


        # given the first random edge, generate a second one where the change would be an actual change
        while True:
            source_2 = choice(edge_source_sample)
            target_2 = choice(list(network.neighbors(source_2)))
            edge_2 = (source_2, target_2)
            if (switch_type == "place" and (len(set(edge_1 + edge_2)) > 2)) or\
               (switch_type != "place" and (len(set(edge_1 + edge_2)) == 4)):
                break
        edge_2 = (source_2, target_2)
        edge_2_attrs = dict(network.edges[edge_2])

        network.remove_edge(*edge_1)
        network.remove_edge(*edge_2)
        if switch_type == "place":
            network.add_edge(*edge_1, **edge_2_attrs)
            network.add_edge(*edge_2, **edge_1_attrs)
        if switch_type == "target":
            network.add_edge(edge_1[0], edge_2[1], **edge_1_attrs)
            network.add_edge(edge_2[0], edge_1[1], **edge_2_attrs)
        if switch_type == "source":
            network.add_edge(edge_2[0], edge_1[1], **edge_1_attrs)
            network.add_edge(edge_1[0], edge_2[1], **edge_2_attrs)

        switches -= 1


def simple_generate_rank_equivalent_network(network: nx.Graph, edge_switch_factor):
    """
    has only one type of switch: target switch. Since each edge has equaly chances of being selected in one direction
    as in the other, this should be the same as tossing a coin to decide which kind of switch to perform.
    :param network:
    :param edge_switch_factor:
    :return:
    """
    switches = int(edge_switch_factor * len(list(network.edges)))
    degree_weighted_nodes = []
    for node in network.nodes:
        degree_weighted_nodes.extend([node] * network.degree(node))

    while switches:
        # randomly select and edge
        source_1 = choice(degree_weighted_nodes)
        target_1 = choice(list(network.neighbors(source_1)))
        edge_1 = (source_1, target_1)

        # Randomly select a second edge that does not share a node with the first
        while True:
            source_2 = choice(degree_weighted_nodes)
            if source_2 not in edge_1:
                target_2 = choice(list(network.neighbors(source_2)))
                if target_2 not in edge_1:
                    break
        edge_2 = (source_2, target_2)

        if (source_1, target_2) in network.edges or (source_2, target_1) in network.edges:
            continue

        network.add_edge(source_1, target_2, **network.edges[edge_1])
        network.add_edge(source_2, target_1, **network.edges[edge_2])
        network.remove_edge(*edge_1)
        network.remove_edge(*edge_2)
        switches -= 1


def randomizing_worker_main(original_network_file, network_loader_class, queue):
    network_loader = network_loader_class(original_network_file)

    while True:
        randomization_request = queue.get(block=True)

        if randomization_request == "HALT":
            break
        output_path = randomization_request["output_path"]
        switch_factor = randomization_request["edge_switch_factor"]
        print(f"now handling request for {output_path}")
        network = network_loader.load_network()
        simple_generate_rank_equivalent_network(network, switch_factor)
        network_loader_class.record_network(network, output_path)


def randomize_network(number_of_randomizations, edge_switch_factor, original_network_file,
                      network_loader_class, output_dir):
    queue = Queue()
    randomizers = Pool(12, randomizing_worker_main, (original_network_file, network_loader_class, queue))
    original_network_name = os_path.basename(original_network_file).split(".")[0]
    for i in range(number_of_randomizations):
        queue.put({"edge_switch_factor": edge_switch_factor,
                   "output_path": os_path.join(output_dir, original_network_name + f"_{110 + i}")})

    for i in range(10):
        queue.put("HALT")

    queue.close()
    queue.join_thread()
    randomizers.close()
    randomizers.join()


if __name__ == "__main__":
    randomize_network(400, 10,
                      r"C:\studies\code\NetProp\data\H_sapiens\H_spaiens_nov_2021.net", network_loaders.HumanCovidHybridNetworkLoader,
                      r"C:\studies\code\NetProp\data\randomized_h_sapiens_with_covid")