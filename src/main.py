import sys
from propagation.functions import propagate_from_config, \
                                  SUPPRESSED_NODES_ORDERING_KEYWORD, NETWORK_ORDERING_KEYWORD, PRIOR_SET_ORDERING_KEYWORD


def main():
    config_path = sys.argv[1]
    propagate_from_config(config_path, ordering={SUPPRESSED_NODES_ORDERING_KEYWORD: 1, NETWORK_ORDERING_KEYWORD: 2})


if __name__ == "__main__":
    main()
