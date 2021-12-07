import sys
from propagation.functions import propagate_from_config

# TODO: find reasonable way to enforce set behavior on pyantics models. Currently cannot support sets of models so everything is a list
# reason for this is that .dict() will convert any member of the set to a dict, but a dict isn't hashable, which then crushes the program



def main():
    config_path = sys.argv[1]
    propagate_from_config(config_path)

if __name__ == "__main__":
    main()
