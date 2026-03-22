#!/usr/bin/env python3
"""Convert WRF output files to WPS intermediate format for metgrid.exe."""

import re
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Optional

import h5netcdf
import numpy as np
import typer
from typing_extensions import Annotated

from wrf_to_int import WPSUtils
from wrf_to_int.WPSUtils import MapProjection, write_slab

######################################################
# Constants

DEFAULT_PRESSURE_LEVELS_HPA = [
    1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
    600, 550, 500, 450, 400, 350, 300, 250, 200, 150,
    100, 70, 50, 30, 20, 10,
]


######################################################
# Helpers


def _wrf_attr(nc, key):
    """Extract a scalar value from a WRF netCDF global attribute (may be a 1-element array)."""
    val = nc.attrs[key]
    if hasattr(val, 'item'):
        val = val.item()
    if isinstance(val, bytes):
        return val.decode()
    return val


def unstagger(data, axis):
    """Average adjacent points along a staggered WRF dimension."""
    slices_lo = [slice(None)] * data.ndim
    slices_hi = [slice(None)] * data.ndim
    slices_lo[axis] = slice(None, -1)
    slices_hi[axis] = slice(1, None)
    return (data[tuple(slices_lo)] + data[tuple(slices_hi)]) / 2.0




def _write_interp_slabs(intfile, field_3d, pressure_3d, target_levels_pa, proj, hdate, map_source,
                        WPSname, units, desc):
    """Interpolate a 3D field to pressure levels and write all slabs."""
    interp = interp_to_pressure_levels(field_3d, pressure_3d, target_levels_pa)
    for k, plevel_pa in enumerate(target_levels_pa):
        write_slab(intfile, interp[k], plevel_pa, proj, WPSname, hdate, units, map_source, desc)


######################################################
# WRF file collection


def _extract_domain(filename):
    """Extract the domain string (e.g., 'd01') from a wrfout filename."""
    match = re.search(r'wrfout_(d\d+)_', filename)
    if match:
        return match.group(1)
    return None


def collect_wrf_files(paths, domain):
    """
    Resolve input paths to a sorted list of wrfout file Paths.

    Accepts a list of file paths, or a single directory. When a directory is
    given, globs for wrfout_* files (with or without .nc extension). If
    multiple domains are found and --domain is not specified, raises an error.
    """
    resolved = []
    for p in paths:
        p = Path(p)
        if p.is_file():
            resolved.append(p)
        elif p.is_dir():
            wrf_files = sorted(p.glob('wrfout_*'))
            if not wrf_files:
                print(f'Error: No wrfout files found in {p}', file=sys.stderr)
                raise typer.Exit(code=1)

            # Check for multiple domains
            domains = set()
            for f in wrf_files:
                d = _extract_domain(f.name)
                if d:
                    domains.add(d)

            if len(domains) > 1:
                if domain is None:
                    print(
                        f'Error: Multiple domains found in {p}: {sorted(domains)}. '
                        f'Use --domain to specify which one to process.',
                        file=sys.stderr,
                    )
                    raise typer.Exit(code=1)
                wrf_files = [f for f in wrf_files if f'wrfout_{domain}_' in f.name]
                if not wrf_files:
                    print(f'Error: No wrfout files for domain {domain} in {p}', file=sys.stderr)
                    raise typer.Exit(code=1)
            elif domain is not None:
                wrf_files = [f for f in wrf_files if f'wrfout_{domain}_' in f.name]
                if not wrf_files:
                    print(f'Error: No wrfout files for domain {domain} in {p}', file=sys.stderr)
                    raise typer.Exit(code=1)

            resolved.extend(wrf_files)
        else:
            print(f'Error: {p} is not a file or directory', file=sys.stderr)
            raise typer.Exit(code=1)

    return sorted(resolved)


######################################################
# Projection extraction


def extract_projection(nc):
    """Extract WPS-compatible projection info from a WRF netCDF file."""
    map_proj = _wrf_attr(nc, 'MAP_PROJ')

    # WPS intermediate file expects dx/dy in km, WRF stores in meters
    if map_proj == 1:  # Lambert Conformal
        return MapProjection(
            projType=WPSUtils.Projections.LC,
            startLat=float(nc['XLAT'][0, 0, 0]),
            startLon=float(nc['XLONG'][0, 0, 0]),
            startI=1.0, startJ=1.0,
            deltaLat=0.0, deltaLon=0.0,
            dx=float(_wrf_attr(nc, 'DX')) / 1000.0,
            dy=float(_wrf_attr(nc, 'DY')) / 1000.0,
            truelat1=float(_wrf_attr(nc, 'TRUELAT1')),
            truelat2=float(_wrf_attr(nc, 'TRUELAT2')),
            xlonc=float(_wrf_attr(nc, 'STAND_LON')),
        )
    elif map_proj == 2:  # Polar Stereographic
        return MapProjection(
            projType=WPSUtils.Projections.PS,
            startLat=float(nc['XLAT'][0, 0, 0]),
            startLon=float(nc['XLONG'][0, 0, 0]),
            startI=1.0, startJ=1.0,
            deltaLat=0.0, deltaLon=0.0,
            dx=float(_wrf_attr(nc, 'DX')) / 1000.0,
            dy=float(_wrf_attr(nc, 'DY')) / 1000.0,
            truelat1=float(_wrf_attr(nc, 'TRUELAT1')),
            xlonc=float(_wrf_attr(nc, 'STAND_LON')),
        )
    elif map_proj == 3:  # Mercator
        return MapProjection(
            projType=WPSUtils.Projections.MERC,
            startLat=float(nc['XLAT'][0, 0, 0]),
            startLon=float(nc['XLONG'][0, 0, 0]),
            startI=1.0, startJ=1.0,
            deltaLat=0.0, deltaLon=0.0,
            dx=float(_wrf_attr(nc, 'DX')) / 1000.0,
            dy=float(_wrf_attr(nc, 'DY')) / 1000.0,
            truelat1=float(_wrf_attr(nc, 'TRUELAT1')),
            xlonc=float(_wrf_attr(nc, 'STAND_LON')),
        )
    elif map_proj == 6:  # Lat-Lon
        xlat = nc['XLAT'][0]
        xlong = nc['XLONG'][0]
        deltalat = float(xlat[1, 0] - xlat[0, 0])
        deltalon = float(xlong[0, 1] - xlong[0, 0])
        return MapProjection(
            projType=WPSUtils.Projections.LATLON,
            startLat=float(xlat[0, 0]),
            startLon=float(xlong[0, 0]),
            startI=1.0, startJ=1.0,
            deltaLat=deltalat, deltaLon=deltalon,
        )
    else:
        msg = f'Unsupported WRF MAP_PROJ: {map_proj}'
        raise ValueError(msg)


######################################################
# Time parsing


def parse_wrf_times(nc):
    """Return list of (time_idx, hdate_string) from WRF Times variable."""
    times = nc['Times'][:]
    result = []
    for t in range(times.shape[0]):
        time_str = b''.join(times[t]).decode()
        # WPS hdate format: YYYY-MM-DD_HH:00:00
        hdate = time_str[:13] + ':00:00'
        result.append((t, hdate))
    return result


######################################################
# Vertical interpolation


def interp_to_pressure_levels(field_3d, pressure_3d, target_levels_pa):
    """
    Interpolate a 3D field from WRF eta levels to target pressure levels.

    Uses linear interpolation in ln(pressure) space. Returns NaN for
    below-ground and above-model-top levels.

    Parameters
    ----------
    field_3d : np.ndarray
        Shape (nz, ny, nx), field values on eta levels.
    pressure_3d : np.ndarray
        Shape (nz, ny, nx), full pressure in Pa on eta levels.
    target_levels_pa : array-like
        Target pressure levels in Pa.

    Returns
    -------
    np.ndarray
        Shape (n_target, ny, nx), interpolated field.
    """
    nz, ny, nx = field_3d.shape
    n_target = len(target_levels_pa)
    result = np.full((n_target, ny, nx), np.nan, dtype=np.float64)

    log_target = np.log(np.array(target_levels_pa, dtype=np.float64))
    log_pressure = np.log(pressure_3d)

    # Reshape to columns for processing
    log_p_2d = log_pressure.reshape(nz, -1)  # (nz, ny*nx)
    field_2d = field_3d.reshape(nz, -1)
    result_2d = result.reshape(n_target, -1)

    for col in range(ny * nx):
        log_p_col = log_p_2d[:, col]
        field_col = field_2d[:, col]

        # Sort by log-pressure ascending for np.interp
        sort_idx = np.argsort(log_p_col)
        result_2d[:, col] = np.interp(
            log_target, log_p_col[sort_idx], field_col[sort_idx],
            right=np.nan,
        )

    return result


######################################################
# Variable derivation functions


def compute_relative_humidity(t, qvapor, pressure):
    """Compute relative humidity (%) using Bolton (1980)."""
    es = 611.2 * np.exp(17.67 * (t - 273.15) / (t - 273.15 + 243.5))
    e = qvapor * pressure / (0.622 + qvapor)
    return np.clip(100.0 * e / es, 0.0, 100.0)


def compute_dewpoint(qvapor, pressure):
    """Compute dewpoint temperature (K) using inverse Bolton formula."""
    e = qvapor * pressure / (0.622 + qvapor)
    with np.errstate(divide='ignore', invalid='ignore'):
        ln_ratio = np.log(e / 611.2)
        td = 273.15 + 243.5 * ln_ratio / (17.67 - ln_ratio)
    return td


def get_soil_layer_names(nc):
    """Get WPS-style soil layer names from WRF soil layer depths (DZS)."""
    dzs = np.asarray(nc['DZS'][0, :], dtype=np.float64)
    depths_cm = dzs * 100.0
    boundaries = [0.0]
    for d in depths_cm:
        boundaries.append(boundaries[-1] + d)

    layer_names = []
    for k in range(len(dzs)):
        top = int(round(boundaries[k]))
        bot = int(round(boundaries[k + 1]))
        layer_names.append((f'SM{top:03d}{bot:03d}', f'ST{top:03d}{bot:03d}', k))

    return layer_names


######################################################
# Main processing


def process_timestep(nc, t_idx, hdate, proj, intfile, target_pressure_levels_pa, map_source):
    """Process one WRF timestep, writing fields variable-by-variable to reduce memory."""

    # --- pressure_3d is needed for all 3D interpolations ---
    pressure_3d = (
        np.asarray(nc['P'][t_idx], dtype=np.float64)
        + np.asarray(nc['PB'][t_idx], dtype=np.float64)
    )

    # --- Group 1: TT, RH, SPECHUMD (share T and QVAPOR) ---
    theta = np.asarray(nc['T'][t_idx], dtype=np.float64) + 300.0
    temperature_3d = theta * (pressure_3d / 100000.0) ** 0.2854
    del theta

    _write_interp_slabs(intfile, temperature_3d, pressure_3d, target_pressure_levels_pa,
                        proj, hdate, map_source, 'TT', 'K', 'Temperature')

    qvapor_3d = np.asarray(nc['QVAPOR'][t_idx], dtype=np.float64)

    rh_3d = compute_relative_humidity(temperature_3d, qvapor_3d, pressure_3d)
    del temperature_3d
    _write_interp_slabs(intfile, rh_3d, pressure_3d, target_pressure_levels_pa,
                        proj, hdate, map_source, 'RH', '%', 'Relative humidity')
    del rh_3d

    spechumd_3d = qvapor_3d / (1.0 + qvapor_3d)
    del qvapor_3d
    _write_interp_slabs(intfile, spechumd_3d, pressure_3d, target_pressure_levels_pa,
                        proj, hdate, map_source, 'SPECHUMD', 'kg kg-1', 'Specific humidity')
    del spechumd_3d

    # --- Group 2: GHT ---
    ght_3d = unstagger(
        (np.asarray(nc['PH'][t_idx], dtype=np.float64)
         + np.asarray(nc['PHB'][t_idx], dtype=np.float64)) / 9.81,
        axis=0,
    )
    _write_interp_slabs(intfile, ght_3d, pressure_3d, target_pressure_levels_pa,
                        proj, hdate, map_source, 'GHT', 'm', 'Geopotential height')
    del ght_3d

    # --- Group 3: UU, VV (wind rotation requires both) ---
    u_grid = unstagger(np.asarray(nc['U'][t_idx], dtype=np.float64), axis=2)
    v_grid = unstagger(np.asarray(nc['V'][t_idx], dtype=np.float64), axis=1)

    if 'COSALPHA' in nc and 'SINALPHA' in nc:
        cosa = np.asarray(nc['COSALPHA'][t_idx], dtype=np.float64)[np.newaxis, :, :]
        sina = np.asarray(nc['SINALPHA'][t_idx], dtype=np.float64)[np.newaxis, :, :]
        u_earth = u_grid * cosa + v_grid * sina
        v_earth = -u_grid * sina + v_grid * cosa
        del cosa, sina
    else:
        u_earth = u_grid
        v_earth = v_grid
    del u_grid, v_grid

    _write_interp_slabs(intfile, u_earth, pressure_3d, target_pressure_levels_pa,
                        proj, hdate, map_source, 'UU', 'm s-1', 'U-component of wind')
    del u_earth
    _write_interp_slabs(intfile, v_earth, pressure_3d, target_pressure_levels_pa,
                        proj, hdate, map_source, 'VV', 'm s-1', 'V-component of wind')
    del v_earth

    del pressure_3d

    # --- Surface fields ---
    psfc = np.asarray(nc['PSFC'][t_idx], dtype=np.float64)
    write_slab(intfile, psfc, 200100.0, proj, 'PSFC', hdate, 'Pa', map_source, 'Surface pressure')

    t2 = np.asarray(nc['T2'][t_idx], dtype=np.float64)
    hgt = np.asarray(nc['HGT'][t_idx], dtype=np.float64)

    # PMSL via hypsometric reduction
    gamma = 0.0065
    t_mean = t2 + gamma * hgt / 2.0
    pmsl = psfc * np.exp(9.81 * hgt / (287.05 * t_mean))
    write_slab(intfile, pmsl, 201300.0, proj, 'PMSL', hdate, 'Pa', map_source, 'Mean sea level pressure')

    write_slab(intfile, np.asarray(nc['TSK'][t_idx], dtype=np.float64),
               200100.0, proj, 'SKINTEMP', hdate, 'K', map_source, 'Skin temperature')

    write_slab(intfile, t2, 200100.0, proj, 'TT', hdate, 'K', map_source, 'Temperature')

    # 10m winds (earth-relative)
    u10 = np.asarray(nc['U10'][t_idx], dtype=np.float64)
    v10 = np.asarray(nc['V10'][t_idx], dtype=np.float64)
    if 'COSALPHA' in nc and 'SINALPHA' in nc:
        cosa = np.asarray(nc['COSALPHA'][t_idx], dtype=np.float64)
        sina = np.asarray(nc['SINALPHA'][t_idx], dtype=np.float64)
        write_slab(intfile, u10 * cosa + v10 * sina, 200100.0, proj, 'UU', hdate, 'm s-1', map_source, 'U-component of wind')
        write_slab(intfile, -u10 * sina + v10 * cosa, 200100.0, proj, 'VV', hdate, 'm s-1', map_source, 'V-component of wind')
    else:
        write_slab(intfile, u10, 200100.0, proj, 'UU', hdate, 'm s-1', map_source, 'U-component of wind')
        write_slab(intfile, v10, 200100.0, proj, 'VV', hdate, 'm s-1', map_source, 'V-component of wind')

    q2 = np.asarray(nc['Q2'][t_idx], dtype=np.float64)
    write_slab(intfile, q2 / (1.0 + q2), 200100.0, proj, 'SPECHUMD', hdate, 'kg kg-1', map_source, 'Specific humidity')
    write_slab(intfile, compute_dewpoint(q2, psfc), 200100.0, proj, 'DEWPT', hdate, 'K', map_source, 'Dewpoint temperature')
    write_slab(intfile, compute_relative_humidity(t2, q2, psfc), 200100.0, proj, 'RH', hdate, '%', map_source, 'Relative humidity')

    # LANDSEA: WRF XLAND 1=land, 2=water -> WPS 1=land, 0=water
    xland = np.asarray(nc['XLAND'][t_idx], dtype=np.float64)
    write_slab(intfile, np.where(xland < 1.5, 1.0, 0.0), 200100.0, proj, 'LANDSEA', hdate, '0/1 Flag', map_source, 'Land-sea mask')

    # Optional surface fields
    if 'SEAICE' in nc:
        write_slab(intfile, np.asarray(nc['SEAICE'][t_idx], dtype=np.float64), 200100.0, proj, 'SEAICE', hdate, 'fraction', map_source, 'Sea ice fraction')

    if 'SST' in nc:
        write_slab(intfile, np.asarray(nc['SST'][t_idx], dtype=np.float64), 200100.0, proj, 'SST', hdate, 'K', map_source, 'Sea surface temperature')

    write_slab(intfile, hgt, 200100.0, proj, 'SOILHGT', hdate, 'm', map_source, 'Terrain height')

    if 'SNOW' in nc:
        write_slab(intfile, np.asarray(nc['SNOW'][t_idx], dtype=np.float64), 200100.0, proj, 'SNOW', hdate, 'kg m-2', map_source, 'Water equivalent snow depth')

    if 'SNOWH' in nc:
        write_slab(intfile, np.asarray(nc['SNOWH'][t_idx], dtype=np.float64), 200100.0, proj, 'SNOWH', hdate, 'm', map_source, 'Physical snow depth')

    # --- Soil fields ---
    if 'DZS' in nc and 'SMOIS' in nc and 'TSLB' in nc:
        soil_layers = get_soil_layer_names(nc)
        for sm_name, st_name, k in soil_layers:
            write_slab(intfile, np.asarray(nc['SMOIS'][t_idx, k], dtype=np.float64), 200100.0, proj, sm_name, hdate, 'm3 m-3', map_source, 'Soil moisture')
            write_slab(intfile, np.asarray(nc['TSLB'][t_idx, k], dtype=np.float64), 200100.0, proj, st_name, hdate, 'K', map_source, 'Soil temperature')


######################################################
# CLI

app = typer.Typer()


@app.command()
def main(
    wrfout_paths: Annotated[List[Path], typer.Argument(
        help='One or more wrfout file paths, or a single directory containing wrfout files',
        exists=True,
        resolve_path=True,
    )],
    start_date: Annotated[datetime, typer.Option(
        '--start-date', '-s',
        help='Starting date-time to convert',
        formats=['%Y-%m-%d', '%Y-%m-%dT%H:%M:%S', '%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H', '%Y-%m-%d_%H'],
    )],
    end_date: Annotated[datetime, typer.Option(
        '--end-date', '-e',
        help='Ending date-time to convert',
        formats=['%Y-%m-%d', '%Y-%m-%dT%H:%M:%S', '%Y-%m-%d %H:%M:%S', '%Y-%m-%d %H', '%Y-%m-%d_%H'],
    )],
    hour_interval: Annotated[int, typer.Option(
        '--hour-interval', '-h',
        help='Interval in hours between records to convert',
    )] = 6,
    domain: Annotated[Optional[str], typer.Option(
        '--domain', '-d',
        help='WRF domain to process (e.g., d01, d02). Required when directory has multiple domains.',
    )] = None,
    pressure_levels: Annotated[Optional[str], typer.Option(
        '--pressure-levels', '-l',
        help='Comma-separated pressure levels in hPa (default: standard 26 levels)',
    )] = None,
    output_prefix: Annotated[str, typer.Option(
        '--prefix', '-p',
        help='Prefix for output intermediate files',
    )] = 'WRF',
    variables: Annotated[Optional[str], typer.Option(
        '--variables', '-v',
        help='Comma-separated list of WPS variable names to process',
    )] = None,
):
    """Convert WRF output files to WPS intermediate format for metgrid.exe."""

    # Parse pressure levels
    if pressure_levels is not None:
        target_levels_hpa = sorted([float(x) for x in pressure_levels.split(',')], reverse=True)
    else:
        target_levels_hpa = sorted(DEFAULT_PRESSURE_LEVELS_HPA, reverse=True)
    target_levels_pa = [x * 100.0 for x in target_levels_hpa]

    # Collect wrfout files
    wrf_files = collect_wrf_files(wrfout_paths, domain)

    # Set up date range
    start = start_date.replace(minute=0, second=0, microsecond=0)
    end = end_date.replace(minute=0, second=0, microsecond=0)
    intv = timedelta(hours=hour_interval)

    print(f'Start date:     {start}')
    print(f'End date:       {end}')
    print(f'Hour interval:  {hour_interval}')
    print(f'Pressure levels ({len(target_levels_hpa)}): {[int(x) if x == int(x) else x for x in target_levels_hpa]}')
    print(f'Input files:    {len(wrf_files)}')

    # Build set of target datetimes
    target_datetimes = set()
    curr = start
    while curr <= end:
        target_datetimes.add(curr)
        curr += intv

    # Process each file
    for wrf_file in wrf_files:
        print(f'\nProcessing {wrf_file.name}')

        with h5netcdf.File(str(wrf_file), 'r') as nc:
            proj = extract_projection(nc)
            times = parse_wrf_times(nc)
            map_source = f'WRF output ({wrf_file.name})'

            for t_idx, hdate in times:
                # Parse hdate to datetime for filtering
                dt = datetime.strptime(hdate[:13], '%Y-%m-%d_%H')  # noqa: DTZ007
                if dt not in target_datetimes:
                    continue

                datestr = hdate[:13]
                print(f'  Writing {output_prefix}:{datestr}')

                intfile = WPSUtils.IntermediateFile(output_prefix, datestr)
                process_timestep(nc, t_idx, hdate, proj, intfile, target_levels_pa, map_source)
                intfile.close()

    print('\nDone.')
