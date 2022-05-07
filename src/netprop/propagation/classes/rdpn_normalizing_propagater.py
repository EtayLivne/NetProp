from pathlib import Path

import pandas as pd

from .normalizing_propagater import NormalizingPropagater
from netprop.networks.randomization import randomize_preloaded_network
from netprop.networks.loaders import NetpropNetwork
from netprop.propagation import prop_from_args
from netprop.models import PropagationResultModel
from netprop.models import ArgsModel

class RDPNNormalizingPropagater(NormalizingPropagater):
    _ARTIFACTS_ROOT_DIR = "artifacts"
    _PROPS_ARTIFACTS_DIR = "props"
    _NETWORKS_ARTIFACTS_DIR = "networks"

    def norm_method(self, suppressed_nodes: list[str], artifacts_path: str, **kwargs):
        try:
            use_existing = kwargs.get("use_existing", False)
            num_randomizations = kwargs.get("num_randomizations", 100)
            edge_switch_factor = kwargs.get("edge_switch_factor", 10)
            self._create_artifacts(artifacts_path, use_existing, num_randomizations, edge_switch_factor, suppressed_nodes)

            rand_prop_results = self._rand_scores_df_dict(artifacts_path)
            self._add_p_value_to_state(rand_prop_results)
        finally:
            self._remove_prop_files(artifacts_path)

    def _create_artifacts(self, artifact_path: str, use_existing: bool,
                          num_randomizations: int, edge_switch_factor: float, suppressed_nodes: list[str]) -> None:
        artifact_path = Path(artifact_path)
        props_root = Path(artifact_path / self._ARTIFACTS_ROOT_DIR / self._PROPS_ARTIFACTS_DIR)
        rand_networks_root = Path(artifact_path / self._ARTIFACTS_ROOT_DIR / self._NETWORKS_ARTIFACTS_DIR)

        # results_network = Path(artifact_root_path / "artifacts/results")
        if not use_existing:
            for root_path in [props_root, rand_networks_root]:
                root_path.mkdir(parents=True, exist_ok=True)
            randomize_preloaded_network(self.network, num_randomizations, edge_switch_factor, str(rand_networks_root))

        networks = list([str(r) for r in rand_networks_root.glob("*")])
        prop_args = [
            ArgsModel(prop_id="rand" + str(i), network_init_args=[nw], prior_set=self.prior_set, suppressed_set=suppressed_nodes) for i, nw in enumerate(networks)
        ]
        prop_from_args(prop_args, NetpropNetwork, str(props_root))

    def _merge_df_list_on_col(self, df_list: list[pd.DataFrame], col: str):
        if not df_list:
            return pd.DataFrame()
        df = df_list[0].copy()
        for other_df in df_list[1:]:
            df = df.merge(other_df, on=col)
        return df

    def _rand_scores_df_dict(self, artifacts_path: str) -> dict[str, pd.DataFrame]:
        artifacts_path = Path(artifacts_path)
        rand_props_root = Path(artifacts_path / self._ARTIFACTS_ROOT_DIR / self._PROPS_ARTIFACTS_DIR)
        rand_prop_results = {l: dict() for l in self._state.columns}
        for res_file in rand_props_root.glob("*.json"):
            parsed = PropagationResultModel.parse_file(res_file)
            for liquid_name in self._state.columns:
                rand_prop_results[liquid_name][parsed.id] = parsed.prop_scores_as_series(by_liquids=liquid_name)

        rand_prop_results = {
            liquid: self._merge_df_list_on_col([data.pop(k) for k in list(data.keys())], "nodes").set_index("nodes")
            for liquid, data in rand_prop_results.items()
        }
        return rand_prop_results

    def _add_p_value_to_state(self, rand_prop_results: dict[str, pd.DataFrame]):
        row_length = len(self._state.columns)
        for liquid_name in self._state.columns:
            rand_df = rand_prop_results[liquid_name]
            rdpn_vals = dict()
            for node in self._state.index:
                true_score = self._state.loc[node][liquid_name]
                rand_scores = rand_df.loc[node]
                # counts which percentage of the random results are larger than the true. That's the p value.
                rdpn_vals[node] = 1 - len(rand_scores[rand_scores < true_score]) / row_length

            self._state[f"rdpn_p_value_{liquid_name}"] = pd.Series(rdpn_vals)

    def _remove_prop_files(self, artifact_path: str) -> None:
        prop_root = Path(artifact_path) / self._ARTIFACTS_ROOT_DIR / self._PROPS_ARTIFACTS_DIR
        if prop_root.exists():
            prop_files = prop_root.glob("*.json")
            for f in prop_files:
                f.unlink()

    def _assert_already_exists(self, artifact_path: str):
        root = Path(artifact_path) / self._ARTIFACTS_ROOT_DIR
        prop = root / self._PROPS_ARTIFACTS_DIR
        networks = root / self._NETWORKS_ARTIFACTS_DIR
        return root.exists() and prop.exists() and networks.exists()
