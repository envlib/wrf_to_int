# wrf_to_int

<p align="center">
    <em>Utility to convert wrf output to WPS intermediate files</em>
</p>

[![build](https://github.com/wrf_to_int/workflows/Build/badge.svg)](https://github.com/wrf_to_int/actions)
[![codecov](https://codecov.io/gh/mullenkamp/wrf_to_int/branch/main/graph/badge.svg)](https://codecov.io/gh/mullenkamp/wrf_to_int)
[![PyPI version](https://badge.fury.io/py/wrf_to_int.svg)](https://badge.fury.io/py/wrf_to_int)

---

**Source Code**: <a href="https://github.com/wrf_to_int" target="_blank">https://github.com/wrf_to_int</a>

---
## Overview
This package provides a command line utility that converts wrfout files to WPS intermediate files specifically for use as input to metgrid.exe for a new subdomain that is in a different coordinate system than the original WRF domain.

WRF stores data on native eta/sigma levels using an Arakawa C-grid with staggered variables. This tool handles:

- **Vertical interpolation** from eta levels to pressure levels (linear interpolation in ln(pressure) space)
- **Unstaggering** of U, V, and geopotential height fields
- **Wind rotation** from grid-relative to earth-relative using COSALPHA/SINALPHA
- **Variable derivation** (temperature from potential temperature, geopotential height from PH+PHB, relative humidity, dewpoint, sea level pressure, specific humidity)
- **All standard WRF projections**: Lambert Conformal, Polar Stereographic, Mercator, and Lat-Lon

## Installation

```bash
pip install wrf_to_int
# or
uv add wrf_to_int
```

## Usage

### Basic usage

Convert wrfout files from a directory, specifying a date range:

```bash
wrf_to_int /path/to/wrfout/files/ -s 2023-02-10 -e 2023-02-10_12
```

### Multiple input files

Pass individual wrfout files directly:

```bash
wrf_to_int wrfout_d02_2023-02-10_00:00:00.nc wrfout_d02_2023-02-11_00:00:00.nc -s 2023-02-10 -e 2023-02-11
```

### Options

```
wrf_to_int <wrfout_paths...> --start-date DATE --end-date DATE [options]

Arguments:
  wrfout_paths              One or more wrfout file paths, or a single directory

Required options:
  -s, --start-date DATE     Starting date-time to convert
  -e, --end-date DATE       Ending date-time to convert

Optional:
  -h, --hour-interval N     Interval in hours between records (default: 6)
  -d, --domain DOMAIN       WRF domain to process (e.g., d01, d02). Required
                            when a directory contains files from multiple domains.
  -l, --pressure-levels L   Comma-separated pressure levels in hPa
                            (default: 1000,975,950,...,30,20,10)
  -p, --prefix PREFIX       Prefix for output files (default: WRF)
  -v, --variables VARS      Comma-separated WPS variable names to process
```

Date formats supported: `YYYY-MM-DD`, `YYYY-MM-DD_HH`, `YYYY-MM-DDTHH:MM:SS`.

### Output

One WPS intermediate file per timestep, named `{PREFIX}:{YYYY-MM-DD_HH}` (e.g., `WRF:2023-02-10_00`).

### Domain auto-detection

When a directory contains wrfout files from multiple domains (e.g., `wrfout_d01_*` and `wrfout_d02_*`), you must specify which domain to process with `--domain`:

```bash
wrf_to_int /path/to/wrfout/ -s 2023-02-10 -e 2023-02-10 -d d02
```

If the directory contains only one domain, it is auto-detected.

### Custom pressure levels

By default, the tool interpolates to 26 standard pressure levels. To specify custom levels:

```bash
wrf_to_int /path/to/wrfout/ -s 2023-02-10 -e 2023-02-10 -l 1000,850,700,500,300,200,100
```

### WPS pipeline

After generating intermediate files, use them with metgrid.exe by setting `fg_name` in your `namelist.wps`:

```
&metgrid
 fg_name = '/path/to/output/WRF'
 ...
/
```

Then run `metgrid.exe` as usual. The intermediate files contain all the fields metgrid needs for the target subdomain (pressure-level and surface meteorological fields, soil fields, land-sea mask, etc.).

## WPS intermediate file output

### Pressure-level fields (one slab per level)
| WPS Field | Description | WRF Source |
|-----------|-------------|------------|
| TT | Temperature (K) | T, P, PB (theta to actual T) |
| UU | U-wind (m/s) | U (unstaggered, earth-relative) |
| VV | V-wind (m/s) | V (unstaggered, earth-relative) |
| GHT | Geopotential height (m) | PH, PHB (unstaggered) |
| RH | Relative humidity (%) | T, QVAPOR, P, PB |
| SPECHUMD | Specific humidity (kg/kg) | QVAPOR |

### Surface fields
| WPS Field | Description | WRF Source |
|-----------|-------------|------------|
| PSFC | Surface pressure (Pa) | PSFC |
| PMSL | Mean sea level pressure (Pa) | PSFC, T2, HGT |
| SKINTEMP | Skin temperature (K) | TSK |
| TT | 2m temperature (K) | T2 |
| UU | 10m U-wind (m/s) | U10 (earth-relative) |
| VV | 10m V-wind (m/s) | V10 (earth-relative) |
| DEWPT | 2m dewpoint (K) | Q2, PSFC |
| RH | 2m relative humidity (%) | T2, Q2, PSFC |
| LANDSEA | Land-sea mask (0/1) | XLAND |
| SEAICE | Sea ice fraction | SEAICE |
| SST | Sea surface temperature (K) | SST |
| SOILHGT | Terrain height (m) | HGT |
| SNOW | Snow water equivalent (kg/m2) | SNOW |
| SNOWH | Physical snow depth (m) | SNOWH |
| SM/ST | Soil moisture/temperature | SMOIS, TSLB (per layer from DZS) |

## Implementation notes

### WPS intermediate file format

This tool was initially based on the WPS intermediate file writing code from [era5_to_int](https://github.com/era5_to_int). During testing with metgrid.exe, several issues were discovered in the shared `WPSUtils.py` module that only manifested when writing non-LATLON projections (Lambert Conformal, Polar Stereographic, Mercator). These issues were invisible in era5_to_int because ERA5 data always uses the LATLON projection.

**Projection codes**: The WPS intermediate file format uses different integer codes for projections than the WPS internal codes. The correct mapping for the intermediate file is:

| Projection | File code | WPS internal code |
|---|---|---|
| Lat-Lon | 0 | PROJ_LATLON |
| Mercator | 1 | PROJ_MERC |
| Lambert Conformal | 3 | PROJ_LC |
| Gaussian | 4 | PROJ_GAUSS |
| Polar Stereographic | 5 | PROJ_PS |
| Cassini | 6 | PROJ_CASSINI |

The original `WPSUtils.py` had LC=1, PS=2, MERC=3 (matching WPS internal codes, not the file format codes). Using the wrong projection code causes metgrid to misread the projection-specific record, since each projection has a different number of fields (e.g., LC has 8 floats, MERC has 6). This resulted in garbage values for truelat1 and the error: `Set true latitude 1 for all projections!`

**dx/dy units**: metgrid reads dx and dy from the intermediate file and multiplies by 1000 (i.e., it expects values in **km**, not meters). WRF stores DX/DY in meters, so the values must be divided by 1000 before writing. Writing meters (e.g., 3000.0) results in metgrid interpreting it as 3000 km.

**earth_radius units**: Same convention as dx/dy. metgrid multiplies the earth_radius value from the intermediate file by 1000, so it expects the value in **km** (6371.229), not meters (6371229.0). ERA5 (era5_to_int) writes meters but this doesn't cause issues because the LATLON projection doesn't use earth_radius for grid positioning. For Lambert Conformal and other projected grids, the incorrect earth radius causes projection math errors, leading to missing values during interpolation.

## Development

### Setup environment

We use [UV](https://docs.astral.sh/uv/) to manage the development environment and production build. 

```bash
uv sync
```

### Run unit tests

You can run all the tests with:

```bash
uv run pytest
```

### Format the code

Execute the following commands to apply linting and check typing:

```bash
uv run ruff check .
uv run black --check --diff .
uv run mypy --install-types --non-interactive wrf_to_int
```

To auto-format:

```bash
uv run black .
uv run ruff check --fix .
```

## License

This project is licensed under the terms of the Apache Software License 2.0.
