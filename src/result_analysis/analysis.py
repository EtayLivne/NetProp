import os
from result_analysis.comparisons import *
import matplotlib.pyplot as plt


def just_the_tip(file_path):
    with open(os.path.abspath(file_path), 'r') as handler:
        rankings = json.load(handler)
    l2_distances = rankings['l2_distance_dict']
    top_l1 = l2_distances['top_100_l1_distances']
    top_l2 = l2_distances['top_100_l2_distances']
    top_10_l1 = top_l1[:10]
    top_10_l2 = top_l2[:10]
    fsdhgbf = 7

def percent_in_direct_cov_interactions(results_path, symobl_map_path):
    with open(os.path.abspath(results_path), 'r') as handler:
        rankings = json.load(handler)
    with open(os.path.abspath(symobl_map_path), 'r') as handler:
        symbols = json.load(handler)
    l2_distances = rankings['l2_distance_dict']
    top_l1 = l2_distances['top_100_l1_distances']
    top_l2 = l2_distances['top_100_l2_distances']
    l1_knocked_genes = set()
    l2_knocked_genes = set()
    for i in range(len(top_l2)):
        for gene in top_l1[i]["knockout_gene_ids"]:
            l1_knocked_genes.add(gene)
        for gene in top_l2[i]["knockout_gene_ids"]:
            l2_knocked_genes.add(gene)

    direct_cov_interaction_genes = set()
    for value in symbols.values():
        direct_cov_interaction_genes.add(value)

    l1_exclusive = {gene for gene in l1_knocked_genes if gene not in direct_cov_interaction_genes}
    l2_exclusive = {gene for gene in l2_knocked_genes if gene not in direct_cov_interaction_genes}
    l1_exclusive_density = [len([e for e in top_l1[:i + 1] if any([gene_id in l1_exclusive for gene_id in e["knockout_gene_ids"]])])
                            for i in range(len(top_l1))]
    l2_exclusive_density = [len([e for e in top_l2[:i + 1] if any([gene_id in l2_exclusive for gene_id in e["knockout_gene_ids"]])])
                            for i in range(len(top_l2))]
    # l1_exclusive_density = [len([e for e in top_l1[:i+1] if e["knockout_gene_ids"][0] in l1_exclusive]) for i in range(len(top_l1))]
    # l2_exclusive_density = [len([e for e in top_l2[:i+1] if e["knockout_gene_ids"][0] in l2_exclusive]) for i in range(len(top_l2))]

    plt.title(os.path.basename(results_path).split('.')[0])
    plt.plot(l1_exclusive_density, label="l1 elements not in immediate interaction with cov2")
    plt.plot(l2_exclusive_density, label="l2 elements not in immediate interaction with cov2")
    plt.legend()
    plt.show()



def combined_percent_in_direct_cov_interactions(results_paths, symobl_map_path):

    with open(os.path.abspath(symobl_map_path), 'r') as handler:
        symbols = json.load(handler)

    for results_path in results_paths:
        with open(os.path.abspath(results_path), 'r') as handler:
            rankings = json.load(handler)
        l2_distances = rankings['l2_distance_dict']
        top_l1 = l2_distances['top_100_l1_distances']
        top_l2 = l2_distances['top_100_l2_distances']
        l1_knocked_genes = set()
        l2_knocked_genes = set()
        for i in range(len(top_l2)):
            for gene in top_l1[i]["knockout_gene_ids"]:
                l1_knocked_genes.add(gene)
            for gene in top_l2[i]["knockout_gene_ids"]:
                l2_knocked_genes.add(gene)

        direct_cov_interaction_genes = set()
        for value in symbols.values():
            direct_cov_interaction_genes.add(value)

        l1_exclusive = {gene for gene in l1_knocked_genes if gene not in direct_cov_interaction_genes}
        l2_exclusive = {gene for gene in l2_knocked_genes if gene not in direct_cov_interaction_genes}
        l1_exclusive_density = [len([e for e in top_l1[:i + 1] if any([gene_id in l1_exclusive for gene_id in e["knockout_gene_ids"]])])
                                for i in range(len(top_l1))]
        l2_exclusive_density = [len([e for e in top_l2[:i + 1] if any([gene_id in l2_exclusive for gene_id in e["knockout_gene_ids"]])])
                                for i in range(len(top_l2))]
        # l1_exclusive_density = [len([e for e in top_l1[:i+1] if e["knockout_gene_ids"][0] in l1_exclusive]) for i in range(len(top_l1))]
        # l2_exclusive_density = [len([e for e in top_l2[:i+1] if e["knockout_gene_ids"][0] in l2_exclusive]) for i in range(len(top_l2))]
        measurment_type = "single" if "single" in results_path else "double"
        plt.plot(l1_exclusive_density, label=f'{measurment_type}: l1 elements not in immediate interaction with cov2')
        plt.plot(l2_exclusive_density, label=f'{measurment_type}: l2 elements not in immediate interaction with cov2')

    plt.title("l1 l2 share of non immediate cov interactions in top scores")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # percent_in_direct_cov_interactions(r'../../results/l1_l2_single_knockout_signed_l2_norm.json',
    #                                    r'C:\studies\code\NetProp\data\gene_symbol_to_entrezgeneid.json')
    combined_percent_in_direct_cov_interactions([
        r'../../results/l1_l2_single_knockout_signed_l2_norm.json',
        r'../../results/l1_l2_double_knockout_signed_l2_norm.json'
        ],  r'C:\studies\code\NetProp\data\gene_symbol_to_entrezgeneid.json')
    #just_the_tip(r'../../results/l1_l2_single_knockout_signed_l2_norm.json')

    # compare_rankings(os.path.abspath(r'../../results/l1_l2_double_knockout_signed_l2_norm.json'))

    # compare_multiple_rankings([
    #     os.path.abspath(r'../../results/l1_l2_double_knockout_signed_l2_norm.json'),
    #     os.path.abspath(r'../../results/l1_l2_single_knockout_signed_l2_norm.json')
    # ])

















    # with open(r'../../results/l1_l2_double_knockout_fixed.json', 'r') as handler:
    #     rankings = json.load(handler)
    #
    # for dict_name, ranking_dict in rankings.items():
    #     for category, ranking_list in ranking_dict.items():
    #         for element in ranking_list:
    #             # knocked_genes = element["knocked_genes"]
    #             # element["knocked_genes_entrezids"] = knocked_genes
    #             # element["knocked_gene_ids"] = [int(name) for name in knocked_genes[1:-1].split(',')]
    #             del(element["knocked_genes"])
    #
    # with open(os.path.abspath(r'../../results/l1_l2_double_knockout_fixed.json'), 'w') as handler:
    #     json.dump(rankings, handler, indent=4)

    #
    # for dict_name, ranking_dict in rankings.items():
    #     for list_name, actual_list in ranking_dict.items():
    #         new_list = []
    #         for element in actual_list:
    #             knockouts = list(element.keys())[0]
    #             element[knockouts].update({"knocked_genes": knockouts})
    #             new_list.append(element[knockouts])
    #         ranking_dict[list_name] = list(new_list)
    #
    #
    # with open(os.path.abspath(r'../../results/l1_l2_double_knockout_fixed.json'), 'w') as handler:
    #     json.dump(rankings, handler, indent=4)



    # for dict_name, ranking_dict in rankings.items():
    #     for category, ranking_list in ranking_dict.items():
    #         for element in ranking_list:
    #             # knocked_genes = element["knocked_genes"]
    #             # element["knocked_genes_entrezids"] = knocked_genes
    #             # element["knocked_gene_ids"] = [int(name) for name in knocked_genes[1:-1].split(',')]
    #             del(element["knocked_genes"])
    #
    # with open(os.path.abspath(r'../../results/l1_l2_double_knockout_fixed.json'), 'w') as handler:
    #     json.dump(rankings, handler, indent=4)