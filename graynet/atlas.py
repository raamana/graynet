from __future__ import annotations

from pathlib import Path

from graynet import config_graynet as cfg
from graynet.domain import AtlasInfo
from graynet.parcellate import roi_labels_centroids
from graynet.utils import check_atlas
from graynet.volumetric import volumetric_roi_info


def format_node_label(node) -> str:
    if isinstance(node, str):
        return node

    try:
        as_float = float(node)
    except (TypeError, ValueError):
        return str(node)

    if as_float.is_integer():
        return str(int(as_float))

    return str(as_float)


def resolve_atlas(base_feature: str, atlas_spec, node_size=None) -> AtlasInfo:
    atlas_spec, atlas_name = check_atlas(atlas_spec)

    if base_feature in cfg.features_cortical:
        uniq_rois, centroids, roi_labels = roi_labels_centroids(atlas_spec, node_size)
        ignore_label = cfg.null_roi_name
    elif base_feature in cfg.features_volumetric:
        uniq_rois, centroids, roi_labels = volumetric_roi_info(atlas_spec)
        ignore_label = cfg.null_roi_index
    else:
        raise ValueError(
            f"Unrecognized type of base_feature: {base_feature}. "
            f"Choose one of {cfg.base_feature_list}"
        )

    node_labels = tuple(format_node_label(roi) for roi in uniq_rois)
    formatted_centroids = {
        format_node_label(roi): tuple(float(coord) for coord in value)
        for roi, value in centroids.items()
    }

    return AtlasInfo(
        atlas_spec=atlas_spec,
        atlas_name=atlas_name,
        roi_values=tuple(uniq_rois),
        node_labels=node_labels,
        centroids=formatted_centroids,
        roi_labels=roi_labels,
        ignore_label=ignore_label,
    )


def atlas_identifier(atlas_spec, atlas_name: str) -> str:
    if isinstance(atlas_spec, str):
        return atlas_spec

    try:
        atlas_path = Path(atlas_spec)
    except TypeError:
        return atlas_name

    return str(atlas_path.resolve()) if atlas_path.exists() else atlas_name
