import argparse
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from graynet import config_graynet as cfg
from graynet.domain import RunConfig
from graynet.exporters import export_csv, export_graphml
from graynet.pipeline import run
from graynet.utils import as_path, check_atlas, check_features, check_subjects

try:
    __version__ = version("graynet")
except PackageNotFoundError:
    __version__ = "2.0.0"


def _add_common_args(parser: argparse.ArgumentParser, allow_multi_feature: bool = False) -> None:
    parser.add_argument("-i", "--input-dir", dest="input_dir", required=True)
    parser.add_argument("-s", "--subjects", dest="subject_ids_path", required=False)
    parser.add_argument("-a", "--atlas", dest="atlas", default=cfg.default_atlas)
    parser.add_argument("-o", "--out-dir", dest="out_dir", default=None)
    parser.add_argument(
        "-f",
        "--feature",
        dest="features",
        nargs="+" if allow_multi_feature else 1,
        required=True,
        help="Feature name(s) to process.",
    )
    parser.add_argument(
        "-p",
        "--smoothing-param",
        dest="smoothing_param",
        type=int,
        default=cfg.default_smoothing_param,
    )
    parser.add_argument("-n", "--node-size", dest="node_size", type=int, default=cfg.default_node_size)
    parser.add_argument("-c", "--num-procs", dest="num_procs", type=int, default=cfg.default_num_procs)


def _resolve_subject_ids(subject_ids_path, input_dir: Path, features: tuple[str, ...]) -> tuple[str, ...]:
    if subject_ids_path is not None:
        subject_ids, _, _, _ = check_subjects(subject_ids_path)
        return tuple(str(subject_id) for subject_id in subject_ids)

    for feature in features:
        if feature not in cfg.features_freesurfer:
            raise ValueError(
                "Path to subject ID list must be specified when non-Freesurfer features are processed."
            )

    subject_ids = tuple(
        subdir.name for subdir in input_dir.iterdir() if subdir.is_dir() and (subdir / "surf").is_dir()
    )
    if not subject_ids:
        raise ValueError(f"No usable Freesurfer subjects found under {input_dir}")
    return subject_ids


def _resolve_out_dir(input_dir: Path, out_dir) -> Path:
    resolved = as_path(out_dir)
    if resolved is None:
        resolved = input_dir / "graynet_runs"
    resolved.mkdir(exist_ok=True, parents=True)
    return resolved


def _build_run_config(args) -> RunConfig:
    input_dir = Path(args.input_dir).resolve()
    if not input_dir.exists():
        raise IOError(f"Given input directory does not exist: {input_dir}")

    features = tuple(check_features(args.features))
    subject_ids = _resolve_subject_ids(args.subject_ids_path, input_dir, features)
    out_dir = _resolve_out_dir(input_dir, args.out_dir)
    atlas_spec, _ = check_atlas(args.atlas)
    node_size = args.node_size if args.node_size is not None else None

    if args.command == "edges":
        return RunConfig(
            mode="edges",
            input_dir=input_dir,
            out_dir=out_dir,
            subject_ids=subject_ids,
            atlas=atlas_spec,
            smoothing_param=args.smoothing_param,
            node_size=node_size,
            base_features=(features[0],),
            weight_methods=tuple(args.weight_methods),
            num_bins=args.num_bins,
            edge_range=tuple(args.edge_range) if args.edge_range is not None else None,
            return_results=False,
            num_procs=args.num_procs,
        )

    if args.command == "multiedge":
        edge_range_dict = None
        if args.multi_edge_range is not None:
            expected = 2 * len(features)
            if len(args.multi_edge_range) != expected:
                raise ValueError(
                    f"Expected {expected} values for --multi-edge-range, got {len(args.multi_edge_range)}."
                )
            edge_range_dict = {}
            for index, feature in enumerate(features):
                offset = 2 * index
                edge_range_dict[feature] = tuple(args.multi_edge_range[offset : offset + 2])

        return RunConfig(
            mode="multiedge",
            input_dir=input_dir,
            out_dir=out_dir,
            subject_ids=subject_ids,
            atlas=atlas_spec,
            smoothing_param=args.smoothing_param,
            node_size=node_size,
            base_features=features,
            weight_methods=tuple(args.weight_methods),
            num_bins=args.num_bins,
            edge_range_dict=edge_range_dict if edge_range_dict is not None else cfg.edge_range_predefined,
            summary_stats=tuple(args.summary_stats),
            return_results=False,
            num_procs=args.num_procs,
            overwrite_results=True,
        )

    return RunConfig(
        mode="roi-stats",
        input_dir=input_dir,
        out_dir=out_dir,
        subject_ids=subject_ids,
        atlas=atlas_spec,
        smoothing_param=args.smoothing_param,
        node_size=node_size,
        base_features=(features[0],),
        roi_stats=tuple(args.roi_stats),
        return_results=False,
        num_procs=args.num_procs,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="graynet")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    edges = subparsers.add_parser("edges", help="Compute single-feature edge networks.")
    _add_common_args(edges)
    edges.add_argument("-w", "--weight", dest="weight_methods", nargs="+", required=True)
    edges.add_argument("-b", "--num-bins", dest="num_bins", type=int, default=cfg.default_num_bins)
    edges.add_argument(
        "-e",
        "--edge-range",
        dest="edge_range",
        nargs=2,
        type=float,
        default=cfg.default_edge_range,
        metavar=("MIN", "MAX"),
    )

    multiedge = subparsers.add_parser("multiedge", help="Compute multi-feature edge networks.")
    _add_common_args(multiedge, allow_multi_feature=True)
    multiedge.add_argument("-w", "--weight", dest="weight_methods", nargs="+", required=True)
    multiedge.add_argument("-b", "--num-bins", dest="num_bins", type=int, default=cfg.default_num_bins)
    multiedge.add_argument(
        "-t",
        "--summary-stat",
        dest="summary_stats",
        nargs="+",
        default=list(cfg.multi_edge_summary_func_default),
    )
    multiedge.add_argument(
        "-l",
        "--multi-edge-range",
        dest="multi_edge_range",
        nargs="*",
        type=float,
        default=None,
        help="Pairs of min/max values, one pair per feature.",
    )

    roi_stats = subparsers.add_parser("roi-stats", help="Compute ROI summary statistics.")
    _add_common_args(roi_stats)
    roi_stats.add_argument("-r", "--roi-stat", dest="roi_stats", nargs="+", required=True)

    export = subparsers.add_parser("export", help="Export canonical outputs to other formats.")
    export_sub = export.add_subparsers(dest="export_command", required=True)
    export_graphml_parser = export_sub.add_parser("graphml", help="Export GraphML files from a run.")
    export_graphml_parser.add_argument("--run-dir", dest="run_dir", required=True)
    export_graphml_parser.add_argument("--out-dir", dest="out_dir", default=None)
    export_csv_parser = export_sub.add_parser("csv", help="Export CSV files from a run.")
    export_csv_parser.add_argument("--run-dir", dest="run_dir", required=True)
    export_csv_parser.add_argument("--out-dir", dest="out_dir", default=None)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "export":
        run_dir = Path(args.run_dir).resolve()
        export_dir = export_graphml(run_dir, args.out_dir) if args.export_command == "graphml" else export_csv(
            run_dir, args.out_dir
        )
        print(export_dir)
        return export_dir

    run_dir = run(_build_run_config(args))
    print(run_dir)
    return run_dir
