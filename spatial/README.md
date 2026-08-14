# Spatial data

## `hamilton_wenham_towns.geojson`

Municipal boundary polygons for the towns of **Hamilton** and **Wenham**, Massachusetts (Essex County).

### Source

- **Dataset:** MassGIS Data: Municipalities (TOWNSSURVEY)
- **Feature class used:** `TOWNSSURVEY_POLYM` (multipart polygons; one feature per municipality)
- **Publisher:** MassGIS / Commonwealth of Massachusetts
- **Product page:** https://www.mass.gov/info-details/massgis-data-municipalities
- **Access method:** ArcGIS FeatureServer query (GeoJSON), layer `TOWNSSURVEY_POLYM`  
  https://services.arcgis.com/2gdL2gxYNFY2TOUb/arcgis/rest/services/MassMunicipalities/FeatureServer/4
- **Filter:** `TOWN IN ('HAMILTON','WENHAM')`
- **CRS:** WGS 84 (EPSG:4326), via `outSR=4326`
- **Downloaded:** 2026-08-08

### Notes

MassGIS describes this layer as the most accurate representation of Massachusetts municipal boundaries based on the legislatively approved record of municipal boundaries. Authoritative determination of boundary locations can only be provided by a licensed land surveyor.

MassGIS regularly updates the statewide layer; this file is a snapshot of Hamilton and Wenham as of the download date above.

## `north_shore_towns.geojson`

Municipal boundary polygons for **Hamilton**, **Wenham**, **Topsfield**, **Ipswich**, **Rowley**, and **Essex**, Massachusetts (Essex County).

### Source

Same MassGIS Municipalities / `TOWNSSURVEY_POLYM` source as above.

- **Filter:** `TOWN IN ('HAMILTON','WENHAM','TOPSFIELD','IPSWICH','ROWLEY','ESSEX')`
- **CRS:** WGS 84 (EPSG:4326)
- **Downloaded:** 2026-08-08

## `town_population.json`

Census population counts used for per-capita night-lights charts.

- **Fields:** MassGIS `POP2010` and `POP2020` on the municipalities layer (U.S. Census)
- **Method for annual values:** linear interpolation between 2010 and 2020; same annual change extrapolated for 2021–2025
- **Towns:** Hamilton, Wenham, Topsfield, Ipswich, Rowley, Essex
