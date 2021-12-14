import networkx as nx
from random import choice
import os.path as os_path
from multiprocessing import Queue, Pool, cpu_count

import networks.loaders.single as network_loaders


def generate_rank_equivalent_network(network: nx.Graph, edge_switch_factor):
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


def _randomizing_worker_main(network_loader_class, network_loader_args, network_loader_kwargs, queue):
    network_loader = network_loader_class(*network_loader_args, **network_loader_kwargs)
    network = network_loader.load()
    while True:
        randomization_request = queue.get(block=True)

        if randomization_request == "HALT":
            break
        output_path = randomization_request["output_path"]
        switch_factor = randomization_request["edge_switch_factor"]
        print(f"now handling request for {output_path}")
        generate_rank_equivalent_network(network, switch_factor)
        network_loader_class.record_network(network, output_path)


def randomize_network(number_of_randomizations, edge_switch_factor,
                      network_loader_class, network_loadr_args, network_loader_kwargs, output_dir):
    num_workers = cpu_count() - 2
    queue = Queue()
    randomizers = Pool(num_workers, _randomizing_worker_main,
                       (network_loader_class, network_loadr_args, network_loader_kwargs, queue))
    original_network_name ="h_sapiens"
    for i in range(number_of_randomizations):
        queue.put({"edge_switch_factor": edge_switch_factor,
                   "output_path": os_path.join(output_dir, original_network_name + f"_{i}")})

    for i in range(num_workers):
        queue.put("HALT")

    queue.close()
    queue.join_thread()
    randomizers.close()
    randomizers.join()


if __name__ == "__main__":
    randomize_network(500, 10, network_loaders.CombinedHumanCovidNetworkLoader,
                      [r"D:\H_sapiens\H_sapiens_aug_2020.net", r"C:\studies\thesis\code\NetProp\data\HEK293T_SARS-CoV-2 (1).cx", r"D:\symbol_to_entrezgene_2021.json"], dict(),
                      r"D:\randomized_h_sapiens")