from pathlib import Path
from shutil import rmtree
from tempfile import TemporaryDirectory

from netprop.propagation.classes.propagater import Propagater


norm_prop_registry = dict()


class NormalizingPropagater(Propagater):
    def __init_subclass__(cls, **kwargs):
        norm_prop_registry[cls.__name__] = cls

    def propagate(self, prior_set: list[str]=None, suppressed_nodes: list[str]=None, source_confidence: float=None,
                  method="iterative", normalize: bool=True, norm_kwargs: dict=None) -> None:
        super().propagate(prior_set=prior_set, suppressed_nodes=suppressed_nodes, source_confidence=source_confidence)
        if normalize:
            norm_kwargs = norm_kwargs or dict()
            self.normalize(suppressed_nodes, norm_kwargs)

    def normalize(self, suppressed_nodes: list[str], norm_kwargs: dict) -> None:
        artifact_path = norm_kwargs.pop("artifacts_path", None)
        keep_artifacts = norm_kwargs.pop("keep_artifacts", True)
        use_existing = norm_kwargs.get("use_existing", False)
        if artifact_path:
            if use_existing:
                use_existing = self._assert_already_exists(artifact_path)
            elif Path(artifact_path).exists():
                rmtree(artifact_path, ignore_errors=True)

            try:
                self.norm_method(suppressed_nodes, artifact_path, **norm_kwargs)
            finally:
                if Path(artifact_path).exists() and not keep_artifacts:
                    rmtree(artifact_path, ignore_errors=True)

        else:
            with TemporaryDirectory() as temp_artifact_path:
                self.norm_method(suppressed_nodes, temp_artifact_path, **norm_kwargs)

    def norm_method(self, suppressed_nodes: list[str], artifact_path: str, **kwargs):
        raise NotImplemented

    def _assert_already_exists(self, artifact_path: str) -> bool:
        raise NotImplemented









    # def _prop_from_all_norm(self, suppressed_nodes: list[str], artifact_path):
    #     original_prior_set = list(self.prior_set)
    #     self.prior_set = list(self.network.keys())
    #     for n, data in self.network.nodes(data=True):
    #         if "source_of" not in data:
    #             data["source_of"] = []
    #         data["source_of"].append("prop_from_all_p_value")
    #
    #     self.propagate(suppressed_nodes=suppressed_nodes, norm_method=None)
    #     self.prior_set = original_prior_set
    #     for n, data in self.network.nodes(data=True):
    #         data["source_of"] = data.pop(-1)
    #
    # def _mix_prior_set_norm(self, suppressed_nodes: list[str], artifact_path, use_existing: bool=False,
    #                         num_randomizations: int=100):
    #
    #     deg_map = {i: [] for i in set(self.network.degree())}
    #     for node in self.network:
    #         deg_map[self.network.degree[node]].append(node)
    #
    #     prior_sets = [
    #                     [choice(deg_map[self.network.degree[n]]) for n in self.prior_set]
    #                     for _ in range(num_randomizations)
    #                  ]
    #
    #
    #     props_root = Path(artifact_path / "artifacts/props")
    #     if not use_existing:
    #         props_root.mkdir(parents=True, exist_ok=True)
    #         NetpropNetwork.record_network(self.network, str(props_root / 'original_network.json'))


