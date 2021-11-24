import inspect
import sys
import common.network_loaders as network_loaders
from network_loader import BaseNetworkLoader
from propagation_utils import load_config, network_from_config, propagate, multiprocessed_propagate_from_config
from advanced_propagation_utils import propagate_a_lot

# TODO: find reasonable way to enforce set behavior on pyantics models. Currently cannot support sets of models so everything is a list
# reason for this is that .dict() will convert any member of the set to a dict, but a dict isn't hashable, which then crushes the program


def worker_main(ppi_config, output_dir, propagation_queue):
    _name, _cls = 0, 1
    loader_classes = {item[_name]: item[_cls] for item in inspect.getmembers(network_loaders, inspect.isclass)
                      if issubclass(item[_cls], BaseNetworkLoader)}
    network = network_from_config(ppi_config, loader_classes)

    while True:
        propagation_config = propagation_queue.get(block=True)
        if propagation_config == "HALT":
            break
        print(f"Now propagating {propagation_config.id}")
        propagate(network, propagation_config, output_dir)



def main():
    config_path = sys.argv[1]
    propagate_a_lot(config_path)
    # multiprocessed_propagate_from_config(config_path)
    # config = load_config(config_path)
    # propagations_queue = Queue()
    # worker_pool = Pool(10, worker_main, (config.ppi_config, config.output_dir_path, propagations_queue))
    # for propagations_config in config.propagations:
    #     propagations_queue.put(propagations_config)
    # for worker in range(10):
    #     propagations_queue.put("HALT")
    #
    # propagations_queue.close()
    # propagations_queue.join_thread()
    # worker_pool.close()
    # worker_pool.join()




if __name__ == "__main__":
    main()












    ################################################### Choose random knockout nodes for matrix conf ###################################################

    # from advanced_propagation_utils import MatrixConfigModel
    # file_path = r"C:\studies\code\NetProp\src\temp\configurations\matrix\knockouts_matrix_conf.json"
    # conf = MatrixConfigModel.parse_file(file_path)
    #
    # network = network_loaders.HSapeinsNetworkLoader("C:\\studies\\code\\NetProp\\data\\H_sapiens\\H_spaiens_nov_2021.net").load_network()
    # random_nodes = sample(list(network.nodes), 40)
    # conf.suppressed_duplicates = {f"{gene}_knockout": {gene} for gene in random_nodes}
    #
    # with open(file_path, 'w') as json_handler:
    #     json_handler.write(conf.json(indent=4))


    """
    {
  "singular_config_path": "C:\\studies\\code\\NetProp\\src\\temp\\configurations\\simple\\original_h_sapiens_single_covid_node_conf.json",
  "randomized_networks_dir": "C:\\studies\\code\\NetProp\\data\\randomized_h_sapiens_with_covid",
  "suppressed_duplicates": [],
  "output_root_dir": "C:\\studies\\code\\NetProp\\src\\temp\\propagations\\knockouts"
}
    """

    # with open(r"C:\studies\code\NetProp\src\temp\configurations\original_h_sapiens_covid_conf.json", "r") as handler:
    #     conf = json.load(handler)
    # conf["output_dir_path"] = r"C:\studies\code\NetProp\src\temp\propagations\gene_knockouts"
    # propagations = []
    # full_network_propagation = conf["propagations"][0]
    # network = network_loaders.HumanCovidHybridNetworkLoader(conf["ppi_config"]["source_file"]).load_network()
    # for n in [n for n in network.nodes if network.nodes[n]["species_id"] == "human"][1000:4000]:
    #     knockout_propagation = full_network_propagation.copy()
    #     knockout_propagation.update({"suppressed_set": [n], "id": f"{n}_knockout"})
    #     propagations.append(knockout_propagation)
    #
    # conf["propagations"] = propagations
    # print(f"{len(propagations)} single knockouts")
    # with open("temp/configurations/single_knockouts_conf.json", "w") as handler:
    #     json.dump(conf, handler, indent=4)


    # network = network_loaders.HumanCovidHybridNetworkLoader("C:\\studies\\code\\NetProp\\data\\H_sapiens.net").load_network()
    # prior_set = [
    #     {
    #         "id": node,
    #         "source_of": [node]
    #     } for node, data in network.nodes(data=True) if data["species_id"] == "sars-cov-2"
    # ]
    #
    # with open("temp\\conf.json", 'w') as handler:
    #     json.dump({
    #         "ppi_config": {
    #             "id": "mr_dummy",
    #             "source_file": "C:\\studies\\code\\NetProp\\data\\H_sapiens.net",
    #             "loader_class": "HumanCovidHybridNetworkLoader",
    #             "protein_id_class": "entrezgene"
    #         },
    #         "global_propagation_params": {
    #             "prior_set_confidence": 0.8,
    #             "halt_conditions": {
    #                 "min_gap": 1e-4
    #             }
    #         },
    #         "propagations": [
    #             {
    #                 "id": "all_covid_sources",
    #                 "prior_set": prior_set
    #             }
    #         ],
    #         "output_dir": "temp\\propagation_results"
    #     }, handler, indent=4)
