import json
from pathlib import Path

import pyarrow.parquet as pyarrow_parquet

from graynet.run_workflow import cli_run, extract, roiwise_stats_indiv


REPO_DIR = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_DIR / "example_data"
VBM_DIR = EXAMPLE_DIR / "volumetric_CAT12"
SUBJECTS = ["CAM_0002_01"]
SUBJECT_LIST = VBM_DIR / "sub_id_list_test.txt"

NUM_ROI_ATLAS = {
    "cat_aal": 122,
    "cat_lpba40": 56,
    "cat_ibsr": 32,
}


def _links_for(atlas_name: str) -> int:
    num_nodes = NUM_ROI_ATLAS[atlas_name]
    return num_nodes * (num_nodes - 1) // 2


def _run_dir(root: Path) -> Path:
    children = [path for path in root.iterdir() if path.is_dir()]
    assert len(children) == 1
    return children[0]


def _read_metadata(run_dir: Path) -> dict:
    return json.loads((run_dir / "run_metadata.json").read_text())


def test_v2_volumetric_edges_api_supports_all_builtin_atlases(tmp_path):
    for atlas_name in ("cat_aal", "cat_lpba40", "cat_ibsr"):
        atlas_out = tmp_path / atlas_name
        results = extract(
            SUBJECTS,
            VBM_DIR,
            base_feature="spm_cat_gmdensity",
            weight_method_list=["manhattan"],
            atlas=atlas_name,
            smoothing_param=10,
            out_dir=atlas_out,
            return_results=True,
            num_procs=1,
        )

        assert sorted(results) == [("manhattan", "CAM_0002_01")]
        assert results[("manhattan", "CAM_0002_01")].size == _links_for(atlas_name)
        run_dir = _run_dir(atlas_out)
        assert pyarrow_parquet.read_table(run_dir / "edges_raw.parquet").num_rows == _links_for(atlas_name)


def test_v2_volumetric_roi_stats_api_matches_node_count(tmp_path):
    results = roiwise_stats_indiv(
        SUBJECTS,
        VBM_DIR,
        base_feature="spm_cat_gmdensity",
        chosen_roi_stats=["median", "mean"],
        atlas="cat_aal",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert sorted(results) == [("mean", "CAM_0002_01"), ("median", "CAM_0002_01")]
    for roi_values in results.values():
        assert roi_values.size == NUM_ROI_ATLAS["cat_aal"]

    run_dir = _run_dir(tmp_path)
    metadata = _read_metadata(run_dir)
    assert metadata["roi_stat_names"] == ["median", "mean"]


def test_v2_volumetric_cli_supports_builtin_atlas_names_and_metadata(tmp_path):
    run_dir = cli_run(
        [
            "edges",
            "-i",
            str(VBM_DIR),
            "-s",
            str(SUBJECT_LIST),
            "-f",
            "spm_cat_gmdensity",
            "-w",
            "manhattan",
            "-a",
            "cat_aal",
            "-o",
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    metadata = _read_metadata(Path(run_dir))
    assert metadata["atlas_name"] == "cat_aal"
    assert len(metadata["node_labels"]) == NUM_ROI_ATLAS["cat_aal"]
    assert (Path(run_dir) / "edges_raw.parquet").exists()


def test_v2_volumetric_cli_default_export_path_is_graynet_runs(tmp_path):
    run_dir = cli_run(
        [
            "edges",
            "-i",
            str(VBM_DIR),
            "-s",
            str(SUBJECT_LIST),
            "-f",
            "spm_cat_gmdensity",
            "-w",
            "manhattan",
            "-a",
            "cat_aal",
            "-c",
            "1",
        ]
    )

    assert "graynet_runs" in str(run_dir)
    assert Path(run_dir).exists()
