#!/usr/bin/env python3
"""Rebuild street-map overlays so pixels within 300 m of the town line are shown.

The drawn boundary stays the MassGIS town outline. Color stretch (vmin/vmax)
stays the 2nd–98th percentile already stored in meta.json, which was computed
from pixel centers inside the town.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
from matplotlib.path import Path as MplPath
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parent / "data" / "light_NASA"
GEO = ROOT / "docs" / "assets" / "north_shore_towns.geojson"
META_PATH = ROOT / "docs" / "assets" / "town_overlays" / "meta.json"
GRADIENT_PATH = ROOT / "docs" / "assets" / "town_overlays" / "legend_gradient.png"

BUFFER_M = 300.0
PAD_DEG = 0.01
OVERLAY_SCALE = 10
DISPLAY = {"MANCHESTER-BY-THE-SEA": "Manchester"}


def display_name(raw: str) -> str:
    return DISPLAY.get(raw, raw.title())


def polygons(geom: dict) -> list:
    if geom["type"] == "Polygon":
        return [geom["coordinates"]]
    return geom["coordinates"]


def ring_path(ring) -> MplPath:
    return MplPath(np.asarray(ring, dtype=np.float64)[:, :2])


def meters_per_degree(lat_deg: float) -> tuple[float, float]:
    lat = np.deg2rad(lat_deg)
    m_lat = (
        111132.92
        - 559.82 * np.cos(2 * lat)
        + 1.175 * np.cos(4 * lat)
        - 0.0023 * np.cos(6 * lat)
    )
    m_lon = 111412.84 * np.cos(lat) - 93.5 * np.cos(3 * lat) + 0.118 * np.cos(5 * lat)
    return float(m_lon), float(m_lat)


def crop_index(lat, lon, lat0, lat1, lon0, lon1):
    row = np.where((lat >= lat0) & (lat <= lat1))[0]
    col = np.where((lon >= lon0) & (lon <= lon1))[0]
    return int(row.min()), int(row.max() + 1), int(col.min()), int(col.max() + 1)


def load_lat_lon(path: Path):
    with h5py.File(path, "r") as h:
        base = "HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields"
        lat = h[f"{base}/lat"][:]
        lon = h[f"{base}/lon"][:]
    return lat, lon


def load_rad(path: Path):
    with h5py.File(path, "r") as h:
        base = "HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields"
        return h[f"{base}/AllAngle_Composite_Snow_Free"][:].astype(np.float64)


def inside_mask(polys, pts) -> np.ndarray:
    mask = np.zeros(pts.shape[0], dtype=bool)
    for poly in polys:
        inside = ring_path(poly[0]).contains_points(pts)
        for hole in poly[1:]:
            inside &= ~ring_path(hole).contains_points(pts)
        mask |= inside
    return mask


def min_dist_meters(points_lonlat, rings, lat0: float, chunk: int = 800) -> np.ndarray:
    m_lon, m_lat = meters_per_degree(lat0)
    px = points_lonlat[:, 0] * m_lon
    py = points_lonlat[:, 1] * m_lat
    best = np.full(len(px), np.inf, dtype=np.float64)
    for ring in rings:
        xy = np.asarray(ring, dtype=np.float64)[:, :2]
        if len(xy) < 2:
            continue
        if xy.shape[0] == 1 or not np.allclose(xy[0], xy[-1]):
            xy = np.vstack([xy, xy[0]])
        ax = xy[:-1, 0] * m_lon
        ay = xy[:-1, 1] * m_lat
        bx = xy[1:, 0] * m_lon
        by = xy[1:, 1] * m_lat
        for i in range(0, len(ax), chunk):
            sl = slice(i, i + chunk)
            abx = bx[sl] - ax[sl]
            aby = by[sl] - ay[sl]
            ab2 = abx * abx + aby * aby
            ab2 = np.where(ab2 < 1e-6, 1e-6, ab2)
            apx = px[:, None] - ax[sl]
            apy = py[:, None] - ay[sl]
            t = np.clip((apx * abx + apy * aby) / ab2, 0.0, 1.0)
            dx = apx - t * abx
            dy = apy - t * aby
            best = np.minimum(best, np.min(dx * dx + dy * dy, axis=1))
    return np.sqrt(best)


def colorize(values, vmin, vmax, mask, gradient):
    span = vmax - vmin if vmax > vmin else 1.0
    t = np.clip((values - vmin) / span, 0, 1)
    idx = np.clip((t * (gradient.shape[0] - 1)).astype(np.int32), 0, gradient.shape[0] - 1)
    rgb = gradient[idx]
    alpha = np.rint(48 + 192 * np.power(t, 0.6)).astype(np.uint8)
    out = np.zeros(values.shape + (4,), dtype=np.uint8)
    out[..., :3] = rgb
    out[..., 3] = alpha
    bad = (~mask) | (~np.isfinite(values)) | (values < 0)
    out[bad] = 0
    return out


def main() -> None:
    geo = json.loads(GEO.read_text(encoding="utf-8"))
    by_name = {}
    for feat in geo["features"]:
        by_name[display_name(feat["properties"]["TOWN"])] = polygons(feat["geometry"])

    meta = json.loads(META_PATH.read_text(encoding="utf-8"))
    files = sorted(DATA.glob("VNP46A4.A*.h10v04.*.h5"))
    if len(files) != len(meta["years"]):
        raise SystemExit(f"expected {len(meta['years'])} HDF files, found {len(files)}")
    lat, lon = load_lat_lon(files[0])
    dlat = float(lat[1] - lat[0])
    dlon = float(lon[1] - lon[0])
    gradient = np.array(Image.open(GRADIENT_PATH).convert("RGB"))[12]

    prepared = {}
    for name in meta["towns"]:
        polys = by_name[name]
        xy = np.asarray([p[:2] for poly in polys for ring in poly for p in ring], dtype=np.float64)
        pr0, pr1, pc0, pc1 = crop_index(
            lat,
            lon,
            xy[:, 1].min() - PAD_DEG,
            xy[:, 1].max() + PAD_DEG,
            xy[:, 0].min() - PAD_DEG,
            xy[:, 0].max() + PAD_DEG,
        )
        bounds = {
            "south": float(lat[pr1 - 1] + dlat / 2),
            "west": float(lon[pc0] - dlon / 2),
            "north": float(lat[pr0] - dlat / 2),
            "east": float(lon[pc1 - 1] + dlon / 2),
        }
        old = meta["overlays"][name]["bounds"]
        for key in bounds:
            if abs(bounds[key] - old[key]) > 1e-6:
                raise SystemExit(f"{name} bounds changed at {key}: {bounds[key]} vs {old[key]}")
        sub_lon = lon[pc0:pc1]
        sub_lat = lat[pr0:pr1]
        grid_lon, grid_lat = np.meshgrid(sub_lon, sub_lat)
        pts = np.column_stack([grid_lon.ravel(), grid_lat.ravel()])
        inside = inside_mask(polys, pts)
        shown = inside.copy()
        outside = np.flatnonzero(~inside)
        if len(outside):
            dist = min_dist_meters(pts[outside], [ring for poly in polys for ring in poly], float(xy[:, 1].mean()))
            shown[outside] = dist <= BUFFER_M
        shown = shown.reshape(sub_lat.shape[0], sub_lon.shape[0])
        inside = inside.reshape(shown.shape)
        sample = ROOT / "docs" / meta["overlays"][name]["files"][0]["file"]
        old_img = np.array(Image.open(sample).convert("RGBA"))
        if old_img.shape[0] != shown.shape[0] * OVERLAY_SCALE or old_img.shape[1] != shown.shape[1] * OVERLAY_SCALE:
            raise SystemExit(f"{name} image shape {old_img.shape} does not match crop {shown.shape}")
        prepared[name] = {
            "pr0": pr0,
            "pr1": pr1,
            "pc0": pc0,
            "pc1": pc1,
            "shown": shown,
            "inside": inside,
            "added": int((shown & ~inside).sum()),
            "in_town": int(inside.sum()),
        }
        print(f"{name}: in town {prepared[name]['in_town']}, add {prepared[name]['added']} pixels within {BUFFER_M:.0f} m")

    for fp in files:
        year = int(fp.name.split(".")[1][1:5])
        rad = load_rad(fp)
        for name, spec in prepared.items():
            o = meta["overlays"][name]
            sub = rad[spec["pr0"] : spec["pr1"], spec["pc0"] : spec["pc1"]]
            rgba = colorize(sub, o["vmin"], o["vmax"], spec["shown"], gradient)
            rel = f"assets/town_overlays/{name.lower()}/{year}.png"
            path = ROOT / "docs" / rel
            old = np.array(Image.open(path).convert("RGBA"))
            native = old[::OVERLAY_SCALE, ::OVERLAY_SCALE]
            keep = native[..., 3] > 0
            if year == meta["years"][0]:
                lost = int((keep & ~spec["shown"]).sum())
                both = keep & (rgba[..., 3] > 0)
                if both.any():
                    diff = np.abs(rgba[both, :3].astype(np.int16) - native[both, :3].astype(np.int16)).mean()
                else:
                    diff = 0.0
                print(f"  {name} {year}: lost {lost}, rgb mae on shared pixels {diff:.3f}")
                if lost:
                    raise SystemExit(f"{name} would hide {lost} pixels that were already shown")
            rgba[keep] = native[keep]
            img = Image.fromarray(rgba, "RGBA")
            img = img.resize((img.width * OVERLAY_SCALE, img.height * OVERLAY_SCALE), Image.Resampling.NEAREST)
            img.save(path, optimize=True)
        print("wrote", year)

    for name, spec in prepared.items():
        meta["overlays"][name]["n_pixels"] = int(spec["shown"].sum())
        meta["overlays"][name]["n_pixels_in_town"] = spec["in_town"]
        meta["overlays"][name]["buffer_m"] = BUFFER_M
    meta["pixel_select"] = "center inside the town line, or within 300 m of that line"
    META_PATH.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("updated", META_PATH)


if __name__ == "__main__":
    main()
