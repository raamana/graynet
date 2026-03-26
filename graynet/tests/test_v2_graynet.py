import json
from functools import partial
from pathlib import Path

import numpy as np
import pyarrow.parquet as pyarrow_parquet
import pytest
from scipy.stats import trim_mean

from graynet import export_to_nx, extract, extract_multiedge, get_edge_values, load_run, roiwise_stats_indiv
from graynet.cli import main as cli_main


REPO_DIR = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_DIR / "example_data"
FS_DIR = EXAMPLE_DIR / "freesurfer"
FS_SUBJECT_LIST = FS_DIR / "list_subjects.txt"
SUBJECTS = ["subject12345"]

NUM_ROI_ATLAS = {
    "fsaverage": 68,
    "glasser2016": 360,
}


def _links_for(atlas_name: str) -> int:
    num_nodes = NUM_ROI_ATLAS[atlas_name]
    return num_nodes * (num_nodes - 1) // 2


def _single_run_dir(root: Path) -> Path:
    children = [path for path in root.iterdir() if path.is_dir()]
    assert len(children) == 1
    return children[0]


def _read_metadata(run_dir: Path) -> dict:
    return json.loads((run_dir / "run_metadata.json").read_text())


def test_v2_edges_api_matches_raw_schema_and_metadata(tmp_path):
    results = extract(
        SUBJECTS,
        FS_DIR,
        base_feature="freesurfer_thickness",
        weight_method_list=["manhattan", "euclidean"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert sorted(results) == [
        ("euclidean", "subject12345"),
        ("manhattan", "subject12345"),
    ]
    for edge_vector in results.values():
        assert edge_vector.size == _links_for("fsaverage")

    run_dir = _single_run_dir(tmp_path)
    raw_table = pyarrow_parquet.read_table(run_dir / "edges_raw.parquet")
    metadata = _read_metadata(run_dir)

    assert raw_table.num_rows == 2 * _links_for("fsaverage")
    assert raw_table.column_names == [
        "subject_id",
        "base_feature",
        "weight_method",
        "u",
        "v",
        "weight",
    ]
    assert metadata["mode"] == "edges"
    assert metadata["weight_methods"] == ["manhattan", "euclidean"]
    assert len(metadata["node_labels"]) == NUM_ROI_ATLAS["fsaverage"]


def test_v2_multiedge_api_writes_raw_and_summary_rows(tmp_path):
    results = extract_multiedge(
        SUBJECTS,
        FS_DIR,
        base_feature_list=["freesurfer_thickness", "freesurfer_curv"],
        weight_method_list=["manhattan", "cosine"],
        summary_stats=["median", "prod"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert sorted(results) == [
        ("cosine", "freesurfer_curv", "subject12345"),
        ("cosine", "freesurfer_thickness", "subject12345"),
        ("manhattan", "freesurfer_curv", "subject12345"),
        ("manhattan", "freesurfer_thickness", "subject12345"),
    ]
    for edge_vector in results.values():
        assert edge_vector.size == _links_for("fsaverage")

    run_dir = _single_run_dir(tmp_path)
    raw_table = pyarrow_parquet.read_table(run_dir / "edges_raw.parquet")
    summary_table = pyarrow_parquet.read_table(run_dir / "edges_summary.parquet")
    metadata = _read_metadata(run_dir)

    assert raw_table.num_rows == 4 * _links_for("fsaverage")
    assert summary_table.num_rows == 4 * _links_for("fsaverage")
    assert metadata["summary_stat_names"] == ["median", "prod"]


def test_v2_roi_stats_api_supports_multiple_stats_and_callable(tmp_path):
    trimmed_mean = partial(trim_mean, proportiontocut=0.05)
    results = roiwise_stats_indiv(
        SUBJECTS,
        FS_DIR,
        base_feature="freesurfer_thickness",
        chosen_roi_stats=["median", np.nanmean, trimmed_mean],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert sorted(results) == [
        ("median", "subject12345"),
        ("nanmean", "subject12345"),
        ("trim_mean", "subject12345"),
    ]
    for roi_values in results.values():
        assert roi_values.size == NUM_ROI_ATLAS["fsaverage"]

    run_dir = _single_run_dir(tmp_path)
    roi_table = pyarrow_parquet.read_table(run_dir / "roi_stats.parquet")
    metadata = _read_metadata(run_dir)
    assert roi_table.num_rows == 3 * NUM_ROI_ATLAS["fsaverage"]
    assert metadata["roi_stat_names"] == ["median", "nanmean", "trim_mean"]


def test_v2_edges_cli_auto_discovers_freesurfer_subjects(tmp_path):
    run_dir = cli_main(
        [
            "edges",
            "-i",
            str(FS_DIR),
            "-f",
            "freesurfer_thickness",
            "-w",
            "manhattan",
            "-a",
            "fsaverage",
            "-o",
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    metadata = _read_metadata(Path(run_dir))
    assert metadata["subject_ids"]
    assert "subject12345" in metadata["subject_ids"]
    assert (Path(run_dir) / "edges_raw.parquet").exists()


def test_v2_multiedge_cli_accepts_explicit_ranges_and_summary_stats(tmp_path):
    run_dir = cli_main(
        [
            "multiedge",
            "-i",
            str(FS_DIR),
            "-s",
            str(FS_SUBJECT_LIST),
            "-f",
            "freesurfer_thickness",
            "freesurfer_curv",
            "-w",
            "manhattan",
            "-t",
            "median",
            "prod",
            "-l",
            "0.0",
            "5.0",
            "-0.3",
            "0.3",
            "-a",
            "fsaverage",
            "-o",
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    metadata = _read_metadata(Path(run_dir))
    assert metadata["edge_range_dict"] == {
        "freesurfer_thickness": [0.0, 5.0],
        "freesurfer_curv": [-0.3, 0.3],
    }
    assert (Path(run_dir) / "edges_summary.parquet").exists()


def test_v2_export_cli_creates_graphml_and_csv_from_run_outputs(tmp_path):
    run_dir = cli_main(
        [
            "edges",
            "-i",
            str(FS_DIR),
            "-s",
            str(FS_SUBJECT_LIST),
            "-f",
            "freesurfer_thickness",
            "-w",
            "manhattan",
            "-a",
            "fsaverage",
            "-o",
            str(tmp_path / "runs"),
            "-c",
            "1",
        ]
    )

    graphml_dir = Path(cli_main(["export", "graphml", "--run-dir", str(run_dir)]))
    csv_dir = Path(cli_main(["export", "csv", "--run-dir", str(run_dir)]))

    assert next(graphml_dir.rglob("*.graphml")).exists()
    assert next(csv_dir.glob("*.csv")).exists()


def test_v2_load_run_helpers_support_filtering_and_graph_export(tmp_path):
    run_dir = extract(
        SUBJECTS,
        FS_DIR,
        base_feature="freesurfer_thickness",
        weight_method_list=["manhattan"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=False,
        num_procs=1,
    )

    edge_data, metadata = load_run(run_dir)
    subset = get_edge_values(
        edge_data,
        subject_id="subject12345",
        base_feature="freesurfer_thickness",
        weight_method="manhattan",
    )
    graph = export_to_nx(subset)

    assert metadata["atlas_name"] == "fsaverage"
    assert subset.table.num_rows == _links_for("fsaverage")
    assert graph.number_of_edges() == _links_for("fsaverage")
    assert graph.number_of_nodes() == NUM_ROI_ATLAS["fsaverage"]


def test_v2_edge_data_iter_subjects_yields_subject_scoped_views(tmp_path):
    run_dir = cli_main(
        [
            "edges",
            "-i",
            str(FS_DIR),
            "-f",
            "freesurfer_thickness",
            "-w",
            "manhattan",
            "-a",
            "fsaverage",
            "-o",
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    edge_data, _ = load_run(run_dir)
    seen = []
    for subject_id, subject_edges in edge_data.iter_subjects():
        seen.append(subject_id)
        assert subject_edges.table.num_rows == _links_for("fsaverage")

    assert "subject12345" in seen
    assert seen == list(edge_data.metadata["subject_ids"])


def test_v2_get_edge_values_accepts_multiple_subject_ids_and_preserves_requested_order(tmp_path):
    run_dir = cli_main(
        [
            "edges",
            "-i",
            str(FS_DIR),
            "-f",
            "freesurfer_thickness",
            "-w",
            "manhattan",
            "-a",
            "fsaverage",
            "-o",
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    edge_data, _ = load_run(run_dir)
    requested_ids = list(reversed(edge_data.metadata["subject_ids"]))
    subset = get_edge_values(
        edge_data,
        subject_id=requested_ids,
        base_feature="freesurfer_thickness",
        weight_method="manhattan",
    )
    matrix = subset.to_ndarray(
        requested_ids,
        base_feature="freesurfer_thickness",
        weight_method="manhattan",
    )

    assert matrix.shape == (len(requested_ids), _links_for("fsaverage"))
    assert np.isfinite(matrix).all()


def test_v2_to_ndarray_requires_sequence_for_subject_order(tmp_path):
    run_dir = extract(
        SUBJECTS,
        FS_DIR,
        base_feature="freesurfer_thickness",
        weight_method_list=["manhattan"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=False,
        num_procs=1,
    )

    edge_data, _ = load_run(run_dir)
    with pytest.raises(TypeError):
        edge_data.to_ndarray("subject12345", base_feature="freesurfer_thickness", weight_method="manhattan")


def test_v2_invalid_edge_range_raises_value_error():
    with pytest.raises(ValueError):
        extract(SUBJECTS, FS_DIR, edge_range=-1)

    with pytest.raises(ValueError):
        extract(SUBJECTS, FS_DIR, edge_range=[1])

    with pytest.raises(ValueError):
        extract(SUBJECTS, FS_DIR, edge_range=[2, 1])

    with pytest.raises(ValueError):
        extract(SUBJECTS, FS_DIR, edge_range=(1, np.nan))


def test_v2_invalid_num_bins_raises_value_error():
    with pytest.raises(ValueError):
        extract(SUBJECTS, FS_DIR, num_bins=np.nan)

    with pytest.raises(ValueError):
        extract(SUBJECTS, FS_DIR, num_bins=np.inf)


def test_v2_empty_subject_list_raises_value_error():
    with pytest.raises(ValueError):
        extract([], FS_DIR)


def test_v2_non_freesurfer_cli_requires_subject_list(tmp_path):
    with pytest.raises(ValueError):
        cli_main(
            [
                "edges",
                "-i",
                str(EXAMPLE_DIR / "volumetric_CAT12"),
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


def test_v2_atlas_subdivision_changes_edge_count(tmp_path):
    results = extract(
        SUBJECTS,
        FS_DIR,
        base_feature="freesurfer_thickness",
        weight_method_list=["manhattan"],
        atlas="fsaverage",
        node_size=1000,
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert results[("manhattan", "subject12345")].size == 273 * 272 // 2
