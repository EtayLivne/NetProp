import json
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



def compare_rankings(rankings_file_path):
    with open(rankings_file_path, 'r') as handler:
        rankings = json.load(handler)

    # this simplistic, deterministic (seedless) hash function is just used to give a unique i.d to groups of proteins. collision chance > 0 but << 1.
    def poly_hash(num_list):
        return pow(sum(num_list), 2) + sum(num_list)


    # TODO dynamic key finding
    for dict_name, ranking_dict in rankings.items():
        top_l1_key = next(k for k in ranking_dict.keys() if "l1" in k)
        top_l2_key = next(k for k in ranking_dict.keys() if "l2" in k)
        top_l1 = ranking_dict[top_l1_key]
        top_l2 = ranking_dict[top_l2_key]

        top_l1_genes = [set.union(*[{poly_hash(top_l1[j]["knockout_gene_ids"])} for j in range(i + 1)]) for i in range(len(top_l1))]
        top_l2_genes = [set.union(*[{poly_hash(top_l2[j]["knockout_gene_ids"])} for j in range(i + 1)]) for i in range(len(top_l2))]
        l1_l2_commonality = [len(top_l1_genes[i].intersection(top_l2_genes[i]))/max(len(top_l1_genes[i]), len(top_l2_genes[i]))
                             for i in range(len(top_l1))]
        plt.plot(l1_l2_commonality, label=dict_name)

    plt.title(os.path.basename(rankings_file_path).split('.')[0])
    plt.legend()
    plt.show()

