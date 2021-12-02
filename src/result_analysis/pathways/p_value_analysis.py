from pathlib import Path
import json
from os.path import isabs

def calc_p_value(propagation_scores, reference_propagation_scores):
    sorted_scores = sorted([reference_propagation_scores[network] - propagation_scores[network] for network in propagation_scores])
    return 1 - sorted_scores.index(reference_propagation_scores["original"] - propagation_scores["original"]) / len(sorted_scores)


def get_significant(pathway_propagations, reference_propagation,  p_value_threshold=0.05, flow_diff_threshold=0.1):
    significants = dict()

    for pathway, propagation_scores in pathway_propagations.items():
        p_value = calc_p_value(propagation_scores, reference_propagation[pathway])
        if p_value <= p_value_threshold and abs(reference_propagation[pathway]["original"] - propagation_scores["original"]) > flow_diff_threshold*reference_propagation[pathway]["original"]:
            significants[pathway] = p_value
        # if True:
        #     significants[pathway] = {"p_value": p_value, "total_diff": abs(reference_propagation[pathway]["original"] - propagation_scores["original"])/reference_propagation[pathway]["original"]}

    return significants


def record_significants(pathway_propagations_dir: str, reference_propagation_file: str, output_path: str):
    if not isabs(reference_propagation_file):
        raise ValueError("reference propagation must be given as an absolute path!")

    output_dict = {}
    files = [str(Path(pathway_file)) for pathway_file in Path(pathway_propagations_dir).glob("*.json")]
    if reference_propagation_file in files:
        files.remove(reference_propagation_file)

    with open(reference_propagation_file, 'r') as handler:
        reference_propagation = json.load(handler)

    for file in files:
        with open(file, 'r') as handler:
            output_dict[Path(file).stem.replace("_pathways", "")] = get_significant(json.load(handler),
                                                                                    reference_propagation)
    with open(output_path, 'w') as handler:
        json.dump(output_dict, handler, indent=4)


if __name__ == "__main__":
    record_significants(r"C:\studies\code\NetProp\src\temp\knockout_pathways_enrichment\high_propagation",
                        r"C:\studies\code\NetProp\src\temp\knockout_pathways_enrichment\high_propagation\merged_covid_source_original_pathways.json",
                        r"C:\studies\code\NetProp\src\temp\significants_results\temp.json")