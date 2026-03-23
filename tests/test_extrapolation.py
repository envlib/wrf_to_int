"""Tests for ECMWF below-ground extrapolation (Trenberth et al. 1993)."""

import numpy as np
import pytest

from wrf_to_int.wrf_to_int import (
    extrapolate_t_below_ground,
    extrapolate_ght_below_ground,
    interp_to_pressure_levels,
    _R_D,
    _G,
)


def _make_column(t_bot, ght_bot, psfc, p_bot, n_levels=10):
    """Build a simple 1-point, n-level atmosphere for testing.

    Returns pressure_3d, temperature_3d, ght_3d, psfc_2d all shaped for a
    single grid column (nz, 1, 1) or (1, 1).
    """
    # Pressure decreasing with altitude from p_bot to ~200 hPa
    pressures = np.linspace(p_bot, 20000.0, n_levels)[::-1]  # top-down: low-p first
    # Actually WRF is bottom-up: index 0 = bottom
    pressures = pressures[::-1]  # now index 0 = bottom (highest pressure)

    pressure_3d = pressures.reshape(n_levels, 1, 1)

    # Temperature with standard lapse rate
    dz_values = -(_R_D * t_bot / _G) * np.log(pressures / p_bot)
    temperatures = t_bot - 0.0065 * dz_values
    temperature_3d = temperatures.reshape(n_levels, 1, 1)

    # Geopotential height from hypsometric equation
    ght_values = np.zeros(n_levels)
    ght_values[0] = ght_bot
    for k in range(1, n_levels):
        t_avg = (temperatures[k - 1] + temperatures[k]) / 2.0
        ght_values[k] = ght_values[k - 1] - (_R_D * t_avg / _G) * np.log(pressures[k] / pressures[k - 1])
    ght_3d = ght_values.reshape(n_levels, 1, 1)

    psfc_2d = np.array([[psfc]])

    return pressure_3d, temperature_3d, ght_3d, psfc_2d


class TestKnownValueSeaLevel:
    """Test 1: Sea-level column with known physical expectations."""

    def test_temperature_increases_below_ground(self):
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=288.0, ght_bot=0.0, psfc=101325.0, p_bot=101000.0,
        )
        target = np.array([105000.0])  # 1050 hPa — below ground

        result = interp_to_pressure_levels(temperature_3d, pressure_3d, target)
        extrapolate_t_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)

        t_extrap = float(result[0, 0, 0])
        assert t_extrap > 288.0, f'Expected T > 288 K below ground, got {t_extrap}'
        assert t_extrap != 288.0, 'T should not be constant-extrapolated'
        # Roughly 3-4 K warmer for ~350m below sea level
        assert 290.0 < t_extrap < 296.0, f'T={t_extrap} outside expected range'

    def test_ght_negative_below_sea_level(self):
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=288.0, ght_bot=0.0, psfc=101325.0, p_bot=101000.0,
        )
        target = np.array([105000.0])

        result = interp_to_pressure_levels(ght_3d, pressure_3d, target)
        extrapolate_ght_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)

        ght_extrap = float(result[0, 0, 0])
        assert ght_extrap < 0.0, f'Expected GHT < 0 m below sea level, got {ght_extrap}'
        assert ght_extrap != 0.0, 'GHT should not be constant-extrapolated'


class TestHydrostaticConsistency:
    """Test 2: Extrapolated T and GHT satisfy the hypsometric equation."""

    def test_hypsometric_relation(self):
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=288.0, ght_bot=100.0, psfc=101325.0, p_bot=101000.0,
        )
        p1, p2 = 105000.0, 107000.0  # Two below-ground levels
        target = np.array([p1, p2])

        t_result = interp_to_pressure_levels(temperature_3d, pressure_3d, target)
        extrapolate_t_below_ground(t_result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)

        ght_result = interp_to_pressure_levels(ght_3d, pressure_3d, target)
        extrapolate_ght_below_ground(ght_result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)

        t1 = float(t_result[0, 0, 0])
        t2 = float(t_result[1, 0, 0])
        z1 = float(ght_result[0, 0, 0])
        z2 = float(ght_result[1, 0, 0])

        # Hypsometric equation: dZ = -(R_d * T_avg / g) * ln(p2/p1)
        t_avg = (t1 + t2) / 2.0
        dz_expected = -(_R_D * t_avg / _G) * np.log(p2 / p1)
        dz_actual = z2 - z1

        # Allow 5% tolerance (the ECMWF polynomial is an approximation)
        np.testing.assert_allclose(dz_actual, dz_expected, rtol=0.05,
                                   err_msg=f'Hypsometric mismatch: dz_actual={dz_actual:.2f}, dz_expected={dz_expected:.2f}')


class TestGeocatComparison:
    """Test 3: Compare against geocat-comp reference implementation."""

    @pytest.fixture
    def geocat_funcs(self):
        try:
            from geocat.comp.interpolation import _temp_extrapolate, _geo_height_extrapolate
            return _temp_extrapolate, _geo_height_extrapolate
        except ImportError:
            pytest.skip('geocat-comp not installed')

    def test_temperature_matches_geocat(self, geocat_funcs):
        import xarray as xr
        _temp_extrapolate, _ = geocat_funcs

        t_bot_val, p_sfc_val, ps_val, phi_sfc_val = 288.0, 101000.0, 101325.0, 0.0
        p_target = 105000.0

        # Our implementation
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=t_bot_val, ght_bot=phi_sfc_val / _G, psfc=ps_val, p_bot=p_sfc_val,
        )
        target = np.array([p_target])
        result = interp_to_pressure_levels(temperature_3d, pressure_3d, target)
        extrapolate_t_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)
        our_t = float(result[0, 0, 0])

        # geocat-comp reference
        geocat_t = float(_temp_extrapolate(
            xr.DataArray(t_bot_val), p_target, xr.DataArray(p_sfc_val),
            xr.DataArray(ps_val), xr.DataArray(phi_sfc_val),
        ))

        np.testing.assert_allclose(our_t, geocat_t, atol=1e-4,
                                   err_msg=f'T mismatch: ours={our_t}, geocat={geocat_t}')

    def test_ght_matches_geocat(self, geocat_funcs):
        import xarray as xr
        _, _geo_height_extrapolate = geocat_funcs

        t_bot_val, p_sfc_val, ps_val = 288.0, 101000.0, 101325.0
        ght_bot_val = 0.0
        phi_sfc_val = ght_bot_val * _G
        p_target = 105000.0

        # Our implementation
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=t_bot_val, ght_bot=ght_bot_val, psfc=ps_val, p_bot=p_sfc_val,
        )
        target = np.array([p_target])
        result = interp_to_pressure_levels(ght_3d, pressure_3d, target)
        extrapolate_ght_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)
        our_ght = float(result[0, 0, 0])

        # geocat-comp reference (returns height in geopotential meters)
        geocat_ght = float(_geo_height_extrapolate(
            xr.DataArray(t_bot_val), p_target, xr.DataArray(p_sfc_val),
            xr.DataArray(ps_val), xr.DataArray(phi_sfc_val),
        ))

        np.testing.assert_allclose(our_ght, geocat_ght, atol=1e-4,
                                   err_msg=f'GHT mismatch: ours={our_ght}, geocat={geocat_ght}')


class TestContinuityAtSurface:
    """Test 4: No discontinuity at the ground interface."""

    def test_temperature_continuous_at_surface(self):
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=288.0, ght_bot=100.0, psfc=101325.0, p_bot=101000.0,
        )
        # Target at exactly PSFC
        target = np.array([101325.0])

        result = interp_to_pressure_levels(temperature_3d, pressure_3d, target)
        t_interp = float(result[0, 0, 0])  # constant-extrapolated value

        extrapolate_t_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)
        t_extrap = float(result[0, 0, 0])

        # At surface, extrapolated value should be very close to interpolated
        # (the ECMWF formula evaluates to ~tstar at p=ps, which should be close to t_bot)
        np.testing.assert_allclose(t_extrap, t_interp, rtol=0.01,
                                   err_msg=f'Discontinuity at surface: interp={t_interp}, extrap={t_extrap}')

    def test_ght_continuous_at_surface(self):
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=288.0, ght_bot=100.0, psfc=101325.0, p_bot=101000.0,
        )
        target = np.array([101325.0])

        result = interp_to_pressure_levels(ght_3d, pressure_3d, target)
        ght_interp = float(result[0, 0, 0])

        extrapolate_ght_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)
        ght_extrap = float(result[0, 0, 0])

        np.testing.assert_allclose(ght_extrap, ght_interp, rtol=0.01,
                                   err_msg=f'Discontinuity at surface: interp={ght_interp}, extrap={ght_extrap}')


class TestHighTerrain:
    """Test 5: High-terrain column with corrections active."""

    def test_high_terrain_temperature_reasonable(self):
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=270.0, ght_bot=3000.0, psfc=70000.0, p_bot=69500.0,
        )
        target = np.array([85000.0])  # 850 hPa — well below this mountain

        result = interp_to_pressure_levels(temperature_3d, pressure_3d, target)
        extrapolate_t_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)

        t_extrap = float(result[0, 0, 0])
        assert t_extrap > 270.0, f'Expected T > 270 K below ground, got {t_extrap}'
        # Should not exceed ~310 K (unrealistic for ~1500m below 3000m surface)
        assert t_extrap < 310.0, f'T={t_extrap} unreasonably high for high-terrain extrapolation'

    def test_high_terrain_ght_reasonable(self):
        pressure_3d, temperature_3d, ght_3d, psfc_2d = _make_column(
            t_bot=270.0, ght_bot=3000.0, psfc=70000.0, p_bot=69500.0,
        )
        target = np.array([85000.0])

        result = interp_to_pressure_levels(ght_3d, pressure_3d, target)
        extrapolate_ght_below_ground(result, pressure_3d, psfc_2d, temperature_3d, ght_3d, target)

        ght_extrap = float(result[0, 0, 0])
        assert ght_extrap < 3000.0, f'Expected GHT < 3000 m below mountain, got {ght_extrap}'
        # Should be roughly 1400-1600m for 850 hPa
        assert ght_extrap > 500.0, f'GHT={ght_extrap} unreasonably low'
