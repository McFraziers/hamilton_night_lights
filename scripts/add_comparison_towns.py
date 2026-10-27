#!/usr/bin/env python3
"""Add comparison towns and draw every polygon part on the existing map window."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
from matplotlib.path import Path as MplPath
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parent / "data" / "light_NASA"
SPATIAL = ROOT / "spatial" / "north_shore_towns.geojson"
DOCS_GEO = ROOT / "docs" / "assets" / "north_shore_towns.geojson"
POP_PATH = ROOT / "spatial" / "town_population.json"
COMPARE_PATH = ROOT / "docs" / "assets" / "compare_timeseries.json"
PERCAPITA_PATH = ROOT / "docs" / "assets" / "percapita_timeseries.json"
OVERLAY_PATH = ROOT / "docs" / "assets" / "compare_overlay.json"
META_PATH = ROOT / "docs" / "assets" / "town_overlays" / "meta.json"
GRADIENT_PATH = ROOT / "docs" / "assets" / "town_overlays" / "legend_gradient.png"
NEW_GEO = [ROOT / "spatial" / "_new_towns.geojson", ROOT / "spatial" / "_manchester.geojson"]

LON0, LON1 = -71.06, -70.68
LAT0, LAT1 = 42.46, 42.78
FRAME_SCALE = 6
OVERLAY_SCALE = 10
PAD_DEG = 0.01
YEARS = list(range(2012, 2026))

DISPLAY = {"MANCHESTER-BY-THE-SEA": "Manchester"}
NEW_TOWNS = ["Manchester", "Gloucester", "Salem", "Middleton"]
NEW_COLORS = {
    "Manchester": "#f48fb1",
    "Gloucester": "#4dd0e1",
    "Salem": "#ff7043",
    "Middleton": "#c5e1a5",
}


def display_name(raw: str) -> str:
    return DISPLAY.get(raw, raw.title())


def polygons(geom: dict) -> list:
    if geom["type"] == "Polygon":
        return [geom["coordinates"]]
    return geom["coordinates"]


def ring_path(ring) -> MplPath:
    return MplPath(np.asarray(ring, dtype=np.float64)[:, :2])


def load_grid(path: Path):
    h = h5py.File(path, "r")
    base = "HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields"
    lat = h[f"{base}/lat"][:]
    lon = h[f"{base}/lon"][:]
    rad = h[f"{base}/AllAngle_Composite_Snow_Free"][:].astype(np.float64)
    return lat, lon, rad


def crop_index(lat, lon, lat0, lat1, lon0, lon1):
    row = np.where((lat >= lat0) & (lat <= lat1))[0]
    col = np.where((lon >= lon0) & (lon <= lon1))[0]
    return int(row.min()), int(row.max() + 1), int(col.min()), int(col.max() + 1)


def svg_ring(ring, sub_lon, sub_lat) -> str:
    parts = []
    for i, (lo, la, *_) in enumerate(ring):
        x = (lo - sub_lon[0]) / (sub_lon[1] - sub_lon[0]) * FRAME_SCALE + FRAME_SCALE / 2
        y = (la - sub_lat[0]) / (sub_lat[1] - sub_lat[0]) * FRAME_SCALE + FRAME_SCALE / 2
        parts.append(("M" if i == 0 else "L") + f"{x:.2f},{y:.2f}")
    parts.append("Z")
    return "".join(parts)


def mask_for(polys, lon, lat):
    pts_xy = []
    for poly in polys:
        for ring in poly:
            pts_xy.extend((p[0], p[1]) for p in ring)
    arr = np.asarray(pts_xy, dtype=np.float64)
    r0, r1, c0, c1 = crop_index(lat, lon, arr[:, 1].min() - 0.02, arr[:, 1].max() + 0.02, arr[:, 0].min() - 0.02, arr[:, 0].max() + 0.02)
    sub_lon = lon[c0:c1]
    sub_lat = lat[r0:r1]
    grid_lon, grid_lat = np.meshgrid(sub_lon, sub_lat)
    pts = np.column_stack([grid_lon.ravel(), grid_lat.ravel()])
    mask = np.zeros(pts.shape[0], dtype=bool)
    for poly in polys:
        inside = ring_path(poly[0]).contains_points(pts)
        for hole in poly[1:]:
            inside &= ~ring_path(hole).contains_points(pts)
        mask |= inside
    return mask.reshape(sub_lat.shape[0], sub_lon.shape[0]), r0, c0


def population_by_year(pop2010: float, pop2020: float) -> dict[str, float]:
    step = (pop2020 - pop2010) / 10.0
    return {str(year): pop2010 + (year - 2010) * step for year in YEARS}


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
    base_geo = json.loads(SPATIAL.read_text(encoding="utf-8"))
    have = {f["properties"]["TOWN"] for f in base_geo["features"]}
    for path in NEW_GEO:
        extra = json.loads(path.read_text(encoding="utf-8"))
        for feat in extra["features"]:
            if feat["properties"]["TOWN"] not in have:
                base_geo["features"].append(feat)
                have.add(feat["properties"]["TOWN"])
    text = json.dumps(base_geo)
    SPATIAL.write_text(text, encoding="utf-8")
    DOCS_GEO.write_text(text, encoding="utf-8")

    by_name = {}
    pops = {}
    for feat in base_geo["features"]:
        raw = feat["properties"]["TOWN"]
        name = display_name(raw)
        by_name[name] = polygons(feat["geometry"])
        if "POP2010" in feat["properties"]:
            pops[name] = (float(feat["properties"]["POP2010"]), float(feat["properties"]["POP2020"]))

    files = sorted(DATA.glob("VNP46A4.A*.h10v04.*.h5"))
    lat, lon, _ = load_grid(files[0])
    r0, r1, c0, c1 = crop_index(lat, lon, LAT0, LAT1, LON0, LON1)
    sub_lat = lat[r0:r1]
    sub_lon = lon[c0:c1]
    if (c1 - c0) * FRAME_SCALE != 546 or (r1 - r0) * FRAME_SCALE != 462:
        raise SystemExit(f"crop changed: {(c1 - c0) * FRAME_SCALE}x{(r1 - r0) * FRAME_SCALE}")

    compare = json.loads(COMPARE_PATH.read_text(encoding="utf-8"))
    towns = list(compare["towns"])
    for name in NEW_TOWNS:
        if name not in towns:
            towns.append(name)
    paths = []
    for name in towns:
        d = "".join(svg_ring(ring, sub_lon, sub_lat) for poly in by_name[name] for ring in poly if len(ring) >= 3)
        paths.append({"town": name, "d": d})
    overlay = {
        "viewBox": [0, 0, (c1 - c0) * FRAME_SCALE, (r1 - r0) * FRAME_SCALE],
        "paths": paths,
        "bbox": compare["bbox"],
    }
    OVERLAY_PATH.write_text(json.dumps(overlay), encoding="utf-8")
    essex_parts = paths[towns.index("Essex")]["d"].count("M")
    print("Essex outline parts", essex_parts)

    masks = {name: mask_for(by_name[name], lon, lat) for name in NEW_TOWNS}
    yearly = {name: [] for name in NEW_TOWNS}
    stacks = {name: [] for name in NEW_TOWNS}
    for fp in files:
        year = int(fp.name.split(".")[1][1:5])
        _, _, rad = load_grid(fp)
        for name, (mask, rr, cc) in masks.items():
            sub = rad[rr : rr + mask.shape[0], cc : cc + mask.shape[1]]
            vals = sub[mask]
            vals = vals[np.isfinite(vals) & (vals >= 0)]
            yearly[name].append(
                {
                    "year": year,
                    "mean": float(np.mean(vals)) if len(vals) else None,
                    "sum": float(np.sum(vals)) if len(vals) else None,
                    "n": int(len(vals)),
                }
            )
            stacks[name].append(vals)
        print("stats", year)

    gradient = np.array(Image.open(GRADIENT_PATH).convert("RGB"))[12]
    meta = json.loads(META_PATH.read_text(encoding="utf-8"))
    dlat = float(lat[1] - lat[0])
    dlon = float(lon[1] - lon[0])
    for name in NEW_TOWNS:
        all_vals = np.concatenate(stacks[name]) if stacks[name] else np.array([0.0])
        vmin = float(np.percentile(all_vals, 2))
        vmax = float(np.percentile(all_vals, 98))
        polys = by_name[name]
        xy = np.asarray([p[:2] for poly in polys for ring in poly for p in ring], dtype=np.float64)
        pr0, pr1, pc0, pc1 = crop_index(
            lat, lon, xy[:, 1].min() - PAD_DEG, xy[:, 1].max() + PAD_DEG, xy[:, 0].min() - PAD_DEG, xy[:, 0].max() + PAD_DEG
        )
        folder = ROOT / "docs" / "assets" / "town_overlays" / name.lower()
        folder.mkdir(parents=True, exist_ok=True)
        files_meta = []
        # mask on this crop, reused
        sub_lon_t = lon[pc0:pc1]
        sub_lat_t = lat[pr0:pr1]
        grid_lon, grid_lat = np.meshgrid(sub_lon_t, sub_lat_t)
        pts = np.column_stack([grid_lon.ravel(), grid_lat.ravel()])
        tmask = np.zeros(pts.shape[0], dtype=bool)
        for poly in polys:
            inside = ring_path(poly[0]).contains_points(pts)
            for hole in poly[1:]:
                inside &= ~ring_path(hole).contains_points(pts)
            tmask |= inside
        tmask = tmask.reshape(sub_lat_t.shape[0], sub_lon_t.shape[0])
        for fp in files:
            year = int(fp.name.split(".")[1][1:5])
            _, _, rad = load_grid(fp)
            sub = rad[pr0:pr1, pc0:pc1]
            rgba = colorize(sub, vmin, vmax, tmask, gradient)
            img = Image.fromarray(rgba, "RGBA")
            img = img.resize((img.width * OVERLAY_SCALE, img.height * OVERLAY_SCALE), Image.Resampling.NEAREST)
            rel = f"assets/town_overlays/{name.lower()}/{year}.png"
            img.save(ROOT / "docs" / rel, optimize=True)
            files_meta.append({"year": year, "file": rel})
        meta["overlays"][name] = {
            "bounds": {
                "south": float(lat[pr1 - 1] + dlat / 2),
                "west": float(lon[pc0] - dlon / 2),
                "north": float(lat[pr0] - dlat / 2),
                "east": float(lon[pc1 - 1] + dlon / 2),
            },
            "vmin": vmin,
            "vmax": vmax,
            "units": "nW cm^-2 sr^-1",
            "n_pixels": int(tmask.sum()),
            "files": files_meta,
        }
        print("overlay", name, "pixels", int(tmask.sum()), "vmax", round(vmax, 2))
    meta["towns"] = towns
    META_PATH.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    compare["towns"] = towns
    compare["colors"].update(NEW_COLORS)
    by_year = {name: {item["year"]: item for item in items} for name, items in yearly.items()}
    for row in compare["series"]:
        for name in NEW_TOWNS:
            row[name] = by_year[name][row["year"]]["mean"]
    COMPARE_PATH.write_text(json.dumps(compare, indent=2), encoding="utf-8")

    per = json.loads(PERCAPITA_PATH.read_text(encoding="utf-8"))
    per["towns"] = towns
    per["colors"].update(NEW_COLORS)
    pop_doc = json.loads(POP_PATH.read_text(encoding="utf-8"))
    for name in NEW_TOWNS:
        p2010, p2020 = pops[name]
        pop_doc["census"][name] = {"POP2010": int(p2010), "POP2020": int(p2020)}
        pop_doc["by_year"][name] = population_by_year(p2010, p2020)
        per["population"]["census"][name] = pop_doc["census"][name]
        per["population"]["by_year"][name] = pop_doc["by_year"][name]
    for row in per["series"]:
        for name in NEW_TOWNS:
            item = by_year[name][row["year"]]
            pop = pop_doc["by_year"][name][str(row["year"])]
            row[name] = (item["sum"] / pop) * 1000.0 if item["sum"] is not None else None
    POP_PATH.write_text(json.dumps(pop_doc, indent=2), encoding="utf-8")
    PERCAPITA_PATH.write_text(json.dumps(per, indent=2), encoding="utf-8")
    print("wrote comparison assets")


if __name__ == "__main__":
    main()
