import math

import pytest

from searchbudget.core import bump_observables as obs
from searchbudget.core import catalogue, lee
from searchbudget.core import public_obs_map as models
from searchbudget.core import yield_model as ym


def test_the_catalogue_is_56_spectra():
    assert len(catalogue.canonical_order()) == 56
    assert len(catalogue.sorted_spectra()) == 56


def test_event_selections_sum_to_104():
    assert catalogue.n_channels(catalogue.canonical_order()) == 104


def test_scan_and_floor_cover_the_same_axes():
    assert set(obs.SCAN) == set(obs.FLOOR) - set(obs.OBS_MERGE)


def test_every_model_points_at_a_known_spectrum():
    for spectra in models.PUBLIC_OBS.values():
        for o in spectra:
            assert obs.canon(o) in obs.SCAN


def test_n_s_is_the_log_of_the_window_over_the_resolution():
    assert obs.n_s(100.0, 1000.0, 0.05) == pytest.approx(math.log(10.0) / 0.05)
    assert obs.n_s(1000.0, 100.0, 0.05) == 0.0


def test_the_bar_grows_with_the_trials_factor():
    assert lee.z_local(1.0) == pytest.approx(5.0)
    assert lee.z_local(4000.0) == pytest.approx(math.sqrt(25 + 2 * math.log(4000.0)))
    assert lee.z_local(1e5) > lee.z_local(1e4) > lee.z_local(1e3)


def test_a_local_five_sigma_is_worth_less_globally():
    assert lee.z_global(1.0) == pytest.approx(5.0)
    assert lee.z_global(7710.0) < 3.0
    assert lee.z_global(1e6) is None


def test_the_closed_form_tracks_the_exact_tail_solution():
    for N in (1e3, 1e4, 1e5):
        assert abs(lee.z_local(N) - lee.z_exact(N)) < 0.06


def test_phi_and_its_inverse_agree():
    for z in (-2.0, 0.0, 1.5, 3.0):
        assert lee.phi_inv(lee.Phi(z)) == pytest.approx(z, abs=1e-6)
    assert lee.p1(0.0) == pytest.approx(0.5)
    assert lee.p1(5.0) == pytest.approx(2.8665e-7, rel=1e-3)


def test_merge_segments_joins_overlaps_only():
    assert lee.merge_segments([(1, 3), (2, 5), (7, 9)]) == [(1, 5), (7, 9)]
    assert lee.merge_segments([(7, 9), (1, 3)]) == [(1, 3), (7, 9)]


def test_a_thin_window_fails_the_fittability_gate():
    fat = ym.gate(200.0, 5000.0, 0.05, ym.weight("jj"))
    thin = ym.gate(200.0, 5000.0, 0.05, ym.weight("ee") * 1e-9)
    assert fat[3] is True
    assert thin[3] is False


def test_the_gate_truncates_at_the_one_event_mass():
    r, w = 0.05, ym.weight("gg")
    hi_scan, _n, _ev, _ok = ym.gate(100.0, 1e5, r, w)
    assert hi_scan == pytest.approx(ym.one_event_mass(r, w))
