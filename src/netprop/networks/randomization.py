import sys
import networkx as nx
from random import choice
import os.path as os_path
from multiprocessing import Queue, Pool, cpu_count

from .loaders import single as network_loaders


def generate_rank_equivalent_network(network: nx.Graph, edge_switch_factor: float) -> None:
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

    print(f"switches: {switches}")

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


def _randomizing_worker_main(network, queue):

    while True:
        randomization_request = queue.get(block=True)

        if randomization_request == "HALT":
            break
        output_path = randomization_request["output_path"]
        switch_factor = randomization_request["edge_switch_factor"]
        print(f"now handling request for {output_path}")
        generate_rank_equivalent_network(network, switch_factor)
        network_loaders.NetpropNetwork.record_network(network, output_path)


def randomize_network(number_of_randomizations: int, edge_switch_factor: float,
                      network_loader_class, network_loader_args, network_loader_kwargs, output_dir: str,
                      max_workers=cpu_count()-2):
    network_loader = network_loader_class(*network_loader_args, **network_loader_kwargs)
    network = network_loader.load()
    randomize_preloaded_network(network, number_of_randomizations, edge_switch_factor, output_dir, max_workers=max_workers)


def randomize_preloaded_network(network: nx.Graph, number_of_randomizations: int, edge_switch_factor: float, output_dir,
                                max_workers=cpu_count()-2):
    queue = Queue()
    num_workers = min(max_workers, number_of_randomizations)
    worker_pool = Pool(num_workers, _randomizing_worker_main, (network, queue))
    original_network_name = "h_sapiens"
    for i in range(number_of_randomizations):
        queue.put({"edge_switch_factor": edge_switch_factor,
                   "output_path": os_path.join(output_dir, original_network_name + f"_{i}")})

    for i in range(num_workers):
        queue.put("HALT")

    queue.close()
    queue.join_thread()
    worker_pool.close()
    worker_pool.join()


if __name__ == "__main__":
    num_randomizations = int(sys.argv[1])
    edge_switch_fator = int(sys.argv[2])
    human_ppi_path = sys.argv[3]
    cov_to_human_ppi_path = sys.argv[4]
    translator_path = sys.argv[5]
    output_dir_path = sys.argv[6]
    randomize_network(num_randomizations, edge_switch_fator, network_loaders.CombinedHumanCovidNetworkLoader,
                      [human_ppi_path, cov_to_human_ppi_path, translator_path], dict(),
                      output_dir_path)