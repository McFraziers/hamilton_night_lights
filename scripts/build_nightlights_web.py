#!/usr/bin/env python3
"""Build web assets: clipped annual night-light frames + timeseries JSON."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
from matplotlib.path import Path as MplPath
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/home/shares/ohi/fisheries_explore/data/light_NASA")
GEO = ROOT / "spatial" / "hamilton_wenham_towns.geojson"
OUT = ROOT / "web" / "assets"
FRAMES = OUT / "frames"

LON0, LON1 = -71.02, -70.68
LAT0, LAT1 = 42.50, 42.72
VMAX = 40.0


def radiance_to_rgba(arr: np.ndarray, vmax: float = VMAX) -> np.ndarray:
    x = np.clip(arr / vmax, 0, 1)
    x = np.power(x, 0.55)
    r = (10 + 245 * np.power(x, 0.9)).astype(np.uint8)
    g = (16 + 200 * np.power(x, 1.1)).astype(np.uint8)
    b = (32 + 80 * (1 - x) + 40 * x).astype(np.uint8)
    r = np.where(x < 0.02, 8, r).astype(np.uint8)
    g = np.where(x < 0.02, 14, g).astype(np.uint8)
    b = np.where(x < 0.02, 28, b).astype(np.uint8)
    a = np.full(arr.shape, 255, dtype=np.uint8)
    return np.dstack([r, g, b, a])


def main() -> None:
    FRAMES.mkdir(parents=True, exist_ok=True)
    with open(GEO) as f:
        gj = json.load(f)

    town_polys: dict[str, list] = {}
    for feat in gj["features"]:
        name = feat["properties"]["TOWN"].title()
        coords = feat["geometry"]["coordinates"]
        if feat["geometry"]["type"] == "Polygon":
            rings = [coords[0]]
        else:
            rings = [poly[0] for poly in coords]
        town_polys[name] = rings

    town_paths = {
        k: [MplPath(np.array(r)[:, :2]) for r in rings]
        for k, rings in town_polys.items()
    }

    files = sorted(DATA.glob("VNP46A4.A*.h10v04.*.h5"))
    series = []
    frame_meta = []

    # Lock crop window from the first file so every year uses identical pixels
    with h5py.File(files[0], "r") as h:
        base = "HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields"
        lat_ref = h[f"{base}/lat"][:]
        lon_ref = h[f"{base}/lon"][:]
    row = np.where((lat_ref >= LAT0) & (lat_ref <= LAT1))[0]
    col = np.where((lon_ref >= LON0) & (lon_ref <= LON1))[0]
    r0, r1 = int(row.min()), int(row.max() + 1)
    c0, c1 = int(col.min()), int(col.max() + 1)
    sub_lat = lat_ref[r0:r1]
    sub_lon = lon_ref[c0:c1]
    Lon, Lat = np.meshgrid(sub_lon, sub_lat)
    pts = np.column_stack([Lon.ravel(), Lat.ravel()])
    scale = 8

    for fp in files:
        year = int(fp.name.split(".")[1][1:5])
        with h5py.File(fp, "r") as h:
            base = "HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields"
            rad = h[f"{base}/AllAngle_Composite_Snow_Free"][:].astype(np.float64)

        sub = rad[r0:r1, c0:c1]

        in_any = np.zeros(pts.shape[0], dtype=bool)
        town_stats = {}
        for tname, paths in town_paths.items():
            m = np.zeros(pts.shape[0], dtype=bool)
            for p in paths:
                m |= p.contains_points(pts)
            m2 = m.reshape(sub.shape)
            in_any |= m
            vals = sub[m2]
            vals = vals[np.isfinite(vals) & (vals >= 0)]
            town_stats[tname] = {
                "mean": float(np.mean(vals)) if len(vals) else None,
                "sum": float(np.sum(vals)) if len(vals) else None,
                "median": float(np.median(vals)) if len(vals) else None,
                "n_pixels": int(len(vals)),
            }

        in_any = in_any.reshape(sub.shape)
        vals_all = sub[in_any]
        vals_all = vals_all[np.isfinite(vals_all) & (vals_all >= 0)]
        series.append(
            {
                "year": year,
                "mean": float(np.mean(vals_all)),
                "sum": float(np.sum(vals_all)),
                "median": float(np.median(vals_all)),
                "n_pixels": int(len(vals_all)),
                "towns": town_stats,
            }
        )

        img = Image.fromarray(radiance_to_rgba(sub, VMAX), "RGBA")
        img = img.resize((img.width * scale, img.height * scale), Image.Resampling.NEAREST)
        out_png = FRAMES / f"{year}.png"
        img.save(out_png, optimize=True)
        frame_meta.append(
            {"year": year, "file": f"assets/frames/{year}.png", "w": img.width, "h": img.height}
        )
        print(year, series[-1]["mean"])

    # Fixed SVG overlay paths (town lines stay put while lights change)
    def lonlat_to_xy(lon_v: float, lat_v: float):
        x = (lon_v - sub_lon[0]) / (sub_lon[1] - sub_lon[0]) * scale + scale / 2
        y = (lat_v - sub_lat[0]) / (sub_lat[1] - sub_lat[0]) * scale + scale / 2
        return x, y

    paths_out = []
    for tname, rings in town_polys.items():
        for ring in rings:
            parts = []
            for i, (lo, la) in enumerate(ring):
                x, y = lonlat_to_xy(lo, la)
                parts.append(("M" if i == 0 else "L") + f"{x:.2f},{y:.2f}")
            parts.append("Z")
            paths_out.append({"town": tname, "d": "".join(parts)})

    overlay = {
        "viewBox": [0, 0, (c1 - c0) * scale, (r1 - r0) * scale],
        "paths": paths_out,
        "bbox": {"lon_min": LON0, "lon_max": LON1, "lat_min": LAT0, "lat_max": LAT1},
    }
    with open(OUT / "town_overlay.json", "w") as f:
        json.dump(overlay, f)

    payload = {
        "product": "VNP46A4",
        "version": "2.0 / Collection 002",
        "layer": "AllAngle_Composite_Snow_Free",
        "units": "nW cm^-2 sr^-1",
        "tile": "h10v04",
        "bbox": {"lon_min": LON0, "lon_max": LON1, "lat_min": LAT0, "lat_max": LAT1},
        "display_vmax": VMAX,
        "series": series,
        "frames": frame_meta,
        "source": "NASA Black Marble (LAADS DAAC)",
        "note": "Annual mean snow-free nighttime radiance inside Hamilton and Wenham town boundaries.",
    }
    with open(OUT / "timeseries.json", "w") as f:
        json.dump(payload, f, indent=2)
    print("wrote", OUT / "timeseries.json")
    print("wrote", OUT / "town_overlay.json")


if __name__ == "__main__":
    main()
