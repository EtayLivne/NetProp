from pathlib import Path

from networks.network_loader import MultiNetworkLoader
from models.config_models import NetworksParametersModel



class FolderLoader(MultiNetworkLoader):
    def __init__(self, loader_class, folder_path, loader_class_args=None, loader_class_kwargs=None, metadata=None):
        self.loader_class = loader_class
        self.folder_path = folder_path
        self.loader_class_args = loader_class_args if loader_class_args else []
        self.loader_class_kwargs = loader_class_kwargs if loader_class_kwargs else dict()
        self.metadata = metadata if metadata else dict()

    def load(self, *args, **kwargs):
        return [NetworksParametersModel(id=file.stem,
                                        multi=False,
                                        loader_class=self.loader_class,
                                        init_args=[str(file)]+self.loader_class_args,
                                        init_kwargs=self.loader_class_kwargs,
                                        metadata=self.metadata)
                for file in Path(self.folder_path).glob("*")]