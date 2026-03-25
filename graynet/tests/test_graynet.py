from pathlib import Path

import pyarrow.parquet as pyarrow_parquet

from graynet.multi_edge import extract_multiedge
from graynet.run_workflow import cli_run, extract, roiwise_stats_indiv


REPO_DIR = Path(__file__).resolve().parents[2]
EXAMPLE_DIR = REPO_DIR / "example_data"
FS_DIR = EXAMPLE_DIR / "freesurfer"
SUBJECTS = ["subject12345"]
FS_SUBJECT_LIST = FS_DIR / "list_subjects.txt"

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


def test_edges_api_writes_run_level_parquet_and_returns_vectors(tmp_path):
    results = extract(
        SUBJECTS,
        FS_DIR,
        base_feature="freesurfer_thickness",
        weight_method_list=["manhattan"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert sorted(results) == [("manhattan", "subject12345")]
    assert results[("manhattan", "subject12345")].size == _links_for("fsaverage")

    run_dir = _single_run_dir(tmp_path)
    raw_table = pyarrow_parquet.read_table(run_dir / "edges_raw.parquet")
    assert raw_table.num_rows == _links_for("fsaverage")
    assert set(raw_table.column_names) == {
        "subject_id",
        "base_feature",
        "weight_method",
        "u",
        "v",
        "weight",
    }
    assert (run_dir / "run_metadata.json").exists()


def test_multiedge_api_writes_raw_and_summary_tables(tmp_path):
    results = extract_multiedge(
        SUBJECTS,
        FS_DIR,
        base_feature_list=["freesurfer_thickness", "freesurfer_curv"],
        weight_method_list=["manhattan"],
        summary_stats=["median"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert sorted(results) == [
        ("manhattan", "freesurfer_curv", "subject12345"),
        ("manhattan", "freesurfer_thickness", "subject12345"),
    ]
    for value in results.values():
        assert value.size == _links_for("fsaverage")

    run_dir = _single_run_dir(tmp_path)
    raw_table = pyarrow_parquet.read_table(run_dir / "edges_raw.parquet")
    summary_table = pyarrow_parquet.read_table(run_dir / "edges_summary.parquet")
    assert raw_table.num_rows == 2 * _links_for("fsaverage")
    assert summary_table.num_rows == _links_for("fsaverage")


def test_roi_stats_api_writes_run_level_parquet_and_returns_vector(tmp_path):
    results = roiwise_stats_indiv(
        SUBJECTS,
        FS_DIR,
        base_feature="freesurfer_thickness",
        chosen_roi_stats=["median"],
        atlas="fsaverage",
        smoothing_param=10,
        out_dir=tmp_path,
        return_results=True,
        num_procs=1,
    )

    assert sorted(results) == ["subject12345"]
    assert results["subject12345"].size == NUM_ROI_ATLAS["fsaverage"]

    run_dir = _single_run_dir(tmp_path)
    roi_table = pyarrow_parquet.read_table(run_dir / "roi_stats.parquet")
    assert roi_table.num_rows == NUM_ROI_ATLAS["fsaverage"]


def test_edges_cli_writes_run_dir(tmp_path):
    run_dir = cli_run(
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
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    assert Path(run_dir).exists()
    assert (Path(run_dir) / "edges_raw.parquet").exists()


def test_multiedge_cli_writes_summary_table(tmp_path):
    run_dir = cli_run(
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
            "-a",
            "fsaverage",
            "-o",
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    assert Path(run_dir).exists()
    assert (Path(run_dir) / "edges_summary.parquet").exists()


def test_roi_stats_cli_writes_roi_parquet(tmp_path):
    run_dir = cli_run(
        [
            "roi-stats",
            "-i",
            str(FS_DIR),
            "-s",
            str(FS_SUBJECT_LIST),
            "-f",
            "freesurfer_thickness",
            "-r",
            "median",
            "-a",
            "fsaverage",
            "-o",
            str(tmp_path),
            "-c",
            "1",
        ]
    )

    assert Path(run_dir).exists()
    assert (Path(run_dir) / "roi_stats.parquet").exists()


def test_export_cli_creates_graphml_and_csv_exports(tmp_path):
    run_dir = cli_run(
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

    graphml_dir = cli_run(["export", "graphml", "--run-dir", str(run_dir)])
    csv_dir = cli_run(["export", "csv", "--run-dir", str(run_dir)])

    assert next(Path(graphml_dir).rglob("*.graphml")).exists()
    assert next(Path(csv_dir).glob("*.csv")).exists()
