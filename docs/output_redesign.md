# Output Redesign

## Goal

Make graynet write one canonical, run-level result format that:

- preserves edgewise numerical accuracy
- can reconstruct the same graphs currently written as GraphML
- scales better to thousands of subjects
- minimizes file proliferation
- keeps optional export paths for users who want CSV or GraphML


## Design Principles

1. Canonical outputs should be tabular and columnar.
2. GraphML and CSV should be optional exports, not the primary storage format.
3. Results should be grouped by run-level settings that apply uniformly to all rows.
4. The schema must preserve enough information to reconstruct NetworkX graphs 1:1.
5. Single-edge, multi-edge, and ROI-stat outputs should share the same run metadata model.


## Run-Level Grouping

Outputs should be organized per run configuration.

A run configuration is defined by values that are effectively constant across the produced rows:

- atlas
- node_size
- smoothing_param
- num_bins
- feature scope
- weight scope
- edge range policy
- null/background handling mode
- package version

In practice, different atlases or node sizes should produce different output datasets.
They should not be mixed into the same parquet file.


## Canonical Files

For a single run configuration, write:

- `run_metadata.json`
- `edges_raw.parquet`
- `edges_summary.parquet` if summaries are requested
- `roi_stats.parquet` if ROI statistics are requested

Optional exports:

- `exports/graphml/...`
- `exports/csv/...`


## `run_metadata.json`

This file should contain:

- `schema_version`
- `graynet_version`
- `created_at`
- `run_type`
  - `single_edge`
  - `multi_edge`
  - `roi_stats`
- `input_dir`
- `subjects`
- `atlas`
- `node_size`
- `smoothing_param`
- `num_bins`
- `base_feature_list`
- `weight_method_list`
- `summary_stats`
- `edge_range`
- `num_nodes`
- `node_labels`
- `node_positions`
- `notes`

`node_labels` should be stored in a stable order.
`node_positions` should contain centroid coordinates when available.


## `edges_raw.parquet`

This is the canonical edge-level table.

Each row represents one weighted edge instance.

### Required columns

- `subject_id`
- `atlas`
- `node_size`
- `smoothing_param`
- `num_bins`
- `base_feature`
- `weight_method`
- `u`
- `v`
- `weight`

### Optional columns

- `edge_range_min`
- `edge_range_max`
- `u_index`
- `v_index`
- `run_id`

### Multi-edge interpretation

For multi-edge runs, `edges_raw.parquet` contains one row per:

- subject
- base feature
- weight method
- edge `(u, v)`

That means the current multi-graph concept becomes a filtered view over this table:

- fix `subject_id`
- fix `weight_method`
- take multiple `base_feature` rows per `(u, v)`


## `edges_summary.parquet`

This table stores summary edges derived from the raw edge table.

Each row represents one summarized edge.

### Required columns

- `subject_id`
- `atlas`
- `node_size`
- `smoothing_param`
- `num_bins`
- `weight_method`
- `summary_stat`
- `u`
- `v`
- `weight`

### Optional columns

- `feature_scope`
- `u_index`
- `v_index`
- `run_id`

For multi-edge runs, this is the table corresponding to today’s summarized GraphML files.


## `roi_stats.parquet`

Each row represents one ROI statistic value.

### Required columns

- `subject_id`
- `atlas`
- `node_size`
- `smoothing_param`
- `base_feature`
- `stat_name`
- `roi`
- `value`

### Optional columns

- `roi_index`
- `run_id`


## Reconstruction Rules

The schema must support exact graph reconstruction.

### Single-edge graph reconstruction

Filter `edges_raw.parquet` by:

- one `subject_id`
- one `base_feature`
- one `weight_method`

Then build `nx.Graph()` from `(u, v, weight)`.

### Multi-edge graph reconstruction

Filter `edges_raw.parquet` by:

- one `subject_id`
- one `weight_method`

Then build `nx.MultiGraph()` and attach:

- `weight`
- `base_feature`

as edge attributes.

### Summary graph reconstruction

Filter `edges_summary.parquet` by:

- one `subject_id`
- one `weight_method`
- one `summary_stat`

Then build `nx.Graph()` from `(u, v, weight)`.

### Node attributes

Node positions should come from `run_metadata.json`, not from repeated edge rows.


## Validation Against Current Implementation

Before switching the writer, compare the new outputs against the baseline generated from the current implementation.

### Required checks

1. Single-edge weights:
   - reconstruct graph from parquet
   - compare all `(u, v, weight)` values to current GraphML/CSV outputs

2. ROI statistics:
   - compare `(roi, value)` pairs to current ROI CSV output

3. Multi-edge graph:
   - reconstruct `MultiGraph`
   - compare all parallel edges by `(u, v, base_feature, weight)`

4. Summary graph:
   - reconstruct summary graph
   - compare `(u, v, weight)` values to current summary GraphML output

5. Edge count:
   - ensure `n * (n - 1) / 2` per graph view

6. Node identity:
   - ensure node labels and centroid metadata match current graph outputs


## Loop Reorganization

The desired orchestration shape is:

1. iterate subjects
2. load one feature once
3. compute all requested weights for that loaded feature
4. append rows to raw edge records
5. derive summary rows after raw rows for the subject are complete
6. write run-level parquet outputs

This avoids rereading the same subject-feature for every weight method while keeping memory bounded.

The preferred internal model is:

- one in-memory raw edge record accumulator per run
- one in-memory ROI-stat accumulator per run
- optional periodic flush if dataset size grows large


## Export Utilities

Provide utility functions rather than primary duplicate outputs.

Suggested API:

- `edges_to_graph(df)`
- `edges_to_multigraph(df)`
- `edges_to_csv(df, path)`
- `edges_to_graphml(df, metadata, path)`
- `roi_stats_to_csv(df, path)`

These should operate on the canonical parquet-loaded dataframes.


## Migration Notes

This is intentionally a format reset.

- old output naming conventions do not need to be preserved
- old helper scripts that reconstruct filenames will need to be updated
- documentation should explicitly describe the new canonical output contract


## Open Decisions

1. Whether to store `u` and `v` as labels only, or as both labels and indices.
   Recommendation: store both.

2. Whether to write one parquet file per run or partition by subject.
   Recommendation: one parquet per run initially.

3. Whether to depend on `pyarrow`.
   Recommendation: yes, for canonical parquet support.

4. Whether to keep GraphML export in the default CLI.
   Recommendation: no; make it opt-in.
