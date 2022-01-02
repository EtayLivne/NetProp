from pathlib import Path

def listify(elem):
    if isinstance(elem, list):
        return elem
    return [elem]


# def attach_to_root(path: str):
#     return str(Path(VOLUME_ROOT) / Path(path))
