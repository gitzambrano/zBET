import math
from pathlib import Path
import sys

# Ensure repository root is on sys.path when running pytest directly
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import pandas as pd
import pytest

from zBET import (
    INDUCED_TORQUE_MODEL,
    FX_COLEMAN,
    FY_COLEMAN,
    INFLOW_MODELS,
    K_IND,
    MODELS,
    OUTPUTS,
    PROFILE_DRAG_MODEL,
    BladePitch,
    BladeSolidity,
    Geometry,
    axial_condition,
    coefficients,
    coleman_kx,
    collective_pitch,
    inflow_gradients,
    plot_results,
    profile_drag_coefficients,
    radial_integrals,
    radial_moments,
    resolve_pitch,
    resolve_solidity,
    run_sweep,
    save_results,
    trim_hover,
)


GEOM = Geometry(0.20, 1000.0, 1.4, 5.7, 0.20, 0.016)

# Golden benchmark reference points for exact regression verification
_GOLDEN_BENCHMARK = {
    "0.0:0.0": {
        "unif": {"CT": 0.006737546134920464, "CQ": 0.000849073029503454, "CQi": 0.0004497130295034529, "CQ0": 0.00039936000000000106, "CH": -1.1156440354875841e-20, "CHi": 0.0, "CH0": -1.1156440354875841e-20, "CY": 0.0, "CMy": 0.0, "CMx": -0.0, "lambda": 0.058041132547805324, "lambda_i": 0.058041132547805324},
        "col": {"CT": 0.006737546134920464, "CQ": 0.000849073029503454, "CQi": 0.0004497130295034529, "CQ0": 0.00039936000000000106, "CH": -1.1156440354875841e-20, "CHi": 0.0, "CH0": -1.1156440354875841e-20, "CY": 0.0, "CMy": 0.0, "CMx": -0.0, "lambda": 0.058041132547805324, "lambda_i": 0.058041132547805324}
    },
    "0.2:0.0": {
        "unif": {"CT": 0.014166362049137127, "CQ": 0.0009524068974002862, "CQi": 0.0005302056341694758, "CQ0": 0.0004222012632308104, "CH": 0.0004252295976127065, "CHi": 0.00019091275901725752, "CH0": 0.00023431683859544898, "CY": -0.0, "CMy": 0.0, "CMx": -0.0035689562049137135, "lambda": 0.03488902759818302, "lambda_i": 0.03488902759818302},
        "col": {"CT": 0.014166362049137127, "CQ": 0.0009524068974002862, "CQi": 0.0005302056341694758, "CQ0": 0.0004222012632308104, "CH": 0.0004252295976127065, "CHi": 0.00019091275901725752, "CH0": 0.00023431683859544898, "CY": -7.615365725934588e-05, "CMy": 0.003072472303638786, "CMx": -0.0035689562049137135, "lambda": 0.03488902759818302, "lambda_i": 0.03488902759818302}
    },
    "0.4:-0.02": {
        "unif": {"CT": 0.024212779734286357, "CQ": 0.000800695080314181, "CQi": 0.000312541767255478, "CQ0": 0.00048815331305870305, "CH": 0.0006058926807034027, "CHi": 0.00011200153222657268, "CH0": 0.00049389114847683, "CY": -0.0, "CMy": 0.0, "CMx": -0.008487032338867138, "lambda": 0.010234058134738, "lambda_i": 0.030234058134738},
        "col": {"CT": 0.024212779734286357, "CQ": 0.000800695080314181, "CQi": 0.000312541767255478, "CQ0": 0.00048815331305870305, "CH": 0.0006058926807034027, "CHi": 0.00011200153222657268, "CH0": 0.00049389114847683, "CY": -0.00036926970733901817, "CMy": 0.0030872141969559863, "CMx": -0.008487032338867138, "lambda": 0.010234058134738, "lambda_i": 0.030234058134738}
    },
    "0.6:0.02": {
        "unif": {"CT": 0.021970906909166478, "CQ": 0.0011143103828352052, "CQi": 0.0005252242778435747, "CQ0": 0.0005890861049916305, "CH": 0.0014337583346541675, "CHi": 0.0006285509740278318, "CH0": 0.0008052073606263356, "CY": -0.0, "CMy": 0.0, "CMx": -0.010427805129860843, "lambda": 0.03828892385647123, "lambda_i": 0.018288923856471233},
        "col": {"CT": 0.021970906909166478, "CQ": 0.0011143103828352052, "CQi": 0.0005252242778435747, "CQ0": 0.0005890861049916305, "CH": 0.0014337583346541675, "CHi": 0.0006285509740278318, "CH0": 0.0008052073606263356, "CY": -2.104761020896054e-05, "CMy": 0.0017975172840317669, "CMx": -0.010427805129860843, "lambda": 0.03828892385647123, "lambda_i": 0.018288923856471233}
    },
}


def test_regression_matches_baseline_benchmark():
    """Verifies that the analytical implementation reproduces the exact numerical baseline
    for uniform and longitudinal Coleman inflow.
    """
    theta0 = 0.12

    for key, data in _GOLDEN_BENCHMARK.items():
        mu_str, muz_str = key.split(":")
        mu = float(mu_str)
        mu_z = float(muz_str)

        # Uniform inflow
        new_unif = coefficients(mu, mu_z, theta0, GEOM, "uniform")
        for k in ["CT", "CQ", "CQi", "CQ0", "CH", "CHi", "CH0", "CY", "CMy", "CMx", "lambda", "lambda_i"]:
            assert new_unif[k] == pytest.approx(data["unif"][k], abs=1e-13, rel=1e-12)

        # Longitudinal Coleman inflow (fy=0)
        new_col = coefficients(mu, mu_z, theta0, GEOM, "coleman_feingold", fy=0.0)
        for k in ["CT", "CQ", "CQi", "CQ0", "CH", "CHi", "CH0", "CY", "CMy", "CMx", "lambda", "lambda_i"]:
            assert new_col[k] == pytest.approx(data["col"][k], abs=1e-13, rel=1e-12)


def test_solidity_three_input_modes_and_equivalence():
    """Verifies that the 3 solidity options (sigma_ref, sigma_geom, and chords)
    correctly compute all 3 solidity metrics and yield identical coefficients for rectangular blades.
    """
    x0 = 0.20
    r = 1.4
    nb = 4
    s_ref_target = 0.20
    s_geom_target = (1.0 - x0) * s_ref_target  # 0.16
    s_thrust_target = (1.0 - x0 ** 3) * s_ref_target  # 0.1984
    c_equiv = s_ref_target * math.pi * r / nb

    sol_1 = resolve_solidity("sigma_ref", sigma_ref=s_ref_target, radius=r, root_cutout=x0)
    sol_2 = resolve_solidity("sigma_geom", sigma_geom=s_geom_target, radius=r, root_cutout=x0)
    sol_3 = resolve_solidity("chords", n_blades=nb, chord_root=c_equiv, chord_tip=c_equiv, radius=r, root_cutout=x0)

    for sol in [sol_1, sol_2, sol_3]:
        assert sol.sigma_ref == pytest.approx(s_ref_target, rel=1e-12)
        assert sol.sigma_geom == pytest.approx(s_geom_target, rel=1e-12)
        assert sol.sigma_thrust == pytest.approx(s_thrust_target, rel=1e-12)

    g1 = Geometry(1000.0, r, 5.7, x0, 0.016, sol_1)
    g2 = Geometry(1000.0, r, 5.7, x0, 0.016, sol_2)
    g3 = Geometry(1000.0, r, 5.7, x0, 0.016, sol_3)

    res1 = coefficients(0.3, 0.0, 0.12, g1, "drees")
    res2 = coefficients(0.3, 0.0, 0.12, g2, "drees")
    res3 = coefficients(0.3, 0.0, 0.12, g3, "drees")

    for k in OUTPUTS:
        assert res2[k] == pytest.approx(res1[k], abs=1e-14)
        assert res3[k] == pytest.approx(res1[k], abs=1e-14)


def test_ndarc_coleman_feingold_has_lateral_ky():
    """Verifies that the official NDARC Coleman-Feingold model includes lateral inflow gradient Ky = -2*mu by default."""
    mu = 0.3
    lam = 0.05
    kx, ky = inflow_gradients(mu, lam, "coleman_feingold")
    assert ky == pytest.approx(-2.0 * mu)
    assert kx == pytest.approx(15.0 * math.pi / 32.0 * (mu / (math.sqrt(mu * mu + lam * lam) + abs(lam))))

    # Compare CMx with and without lateral Ky
    res_with_ky = coefficients(mu, 0.0, 0.12, GEOM, "coleman_feingold", fy=1.0)
    res_no_ky = coefficients(mu, 0.0, 0.12, GEOM, "coleman_feingold", fy=0.0)
    assert res_with_ky["CMx"] < res_no_ky["CMx"]


def test_figure_of_merit_and_lift_to_drag_ratio():
    """Validates Figure of Merit (FoM) in hover, CPair, and effective L/D ratio in forward flight."""
    # Nominal hover (CT = 0.022):
    th_nom = collective_pitch(GEOM, ct_hover_target=0.022)
    hover_nom = coefficients(0.0, 0.0, th_nom, GEOM, "uniform")
    assert hover_nom["FoM"] == pytest.approx(0.756, rel=0.05)
    assert hover_nom["L_D_eff"] == pytest.approx(0.0)
    assert hover_nom["CPair"] == pytest.approx(hover_nom["CQ"])

    # Hover with arbitrary theta0:
    hover = coefficients(0.0, 0.0, 0.12, GEOM, "uniform")
    assert hover["FoM"] > 0.40
    assert hover["L_D_eff"] == pytest.approx(0.0)
    assert hover["CPair"] == pytest.approx(hover["CQ"])

    # Forward flight (mu = 0.3):
    fwd = coefficients(0.3, 0.0, 0.12, GEOM, "coleman_feingold")
    assert fwd["CPair"] == pytest.approx(
        fwd["CQ"] + 0.3 * fwd["CH"] - 0.0 * fwd["CT"], rel=1e-12
    )
    expected_ld = 0.3 * fwd["CT"] / fwd["CPair"]
    assert fwd["L_D_eff"] == pytest.approx(expected_ld, rel=1e-12)
    assert fwd["L_D_eff"] > 0.0


def test_profile_drag_vectorial_low_order_matches_classical_factors():
    """Validates the tangential and vectorial closed forms against Johnson low-order factors."""
    geom = Geometry(
        1000.0, 1.4, 5.7, 0.0, 0.016,
        resolve_solidity("sigma_ref", sigma_ref=0.20, root_cutout=0.0),
    )
    mu = 0.2
    sigma = geom.sigma
    cd0 = geom.cd0

    ct_t, ch_t, cq_t = profile_drag_coefficients(
        mu, 0.0, geom, profile_drag_model="analytical_tangential"
    )
    ct_v, ch_v, cq_v = profile_drag_coefficients(
        mu, 0.0, geom, profile_drag_model="analytical_vectorial"
    )

    assert ct_t == pytest.approx(0.0)
    assert ct_v == pytest.approx(0.0)
    assert ch_t == pytest.approx(sigma * cd0 * mu / 4.0)
    assert ch_v == pytest.approx(3.0 * sigma * cd0 * mu / 8.0)
    assert cq_t == pytest.approx(sigma * cd0 / 8.0 * (1.0 + mu * mu))
    assert cq_v == pytest.approx(sigma * cd0 / 8.0 * (1.0 + 1.5 * mu * mu))

    cp0_vectorial = cq_v + mu * ch_v
    assert cp0_vectorial == pytest.approx(
        sigma * cd0 / 8.0 * (1.0 + 4.5 * mu * mu)
    )


def test_numerical_vectorial_computes_ct0_and_air_power_identity():
    """Checks normal profile drag and the complete translational-work bookkeeping."""
    mu, mu_z = 0.3, 0.03
    result = coefficients(
        mu, mu_z, 0.12, GEOM, "uniform",
        profile_drag_model="numerical_vectorial",
        induced_torque_model="energy_balance",
    )
    assert result["CT0"] < 0.0
    assert result["CPair"] == pytest.approx(
        result["CQ"] + mu * result["CH"] - mu_z * result["CT"], rel=1e-12
    )

    descending = coefficients(
        mu, -mu_z, 0.12, GEOM, "uniform",
        profile_drag_model="numerical_vectorial",
        induced_torque_model="energy_balance",
    )
    assert descending["CT0"] > 0.0


def test_axial_input_sign_conventions_and_collective_modes():
    """Validates axial flow conversions and collective pitch helpers."""
    mu_z, _ = axial_condition(0.2, 10.0, "alpha", GEOM)
    assert mu_z == pytest.approx(-0.2 * np.tan(np.deg2rad(10.0)))
    assert axial_condition(0.2, 0.03, "mu_z", GEOM)[0] == pytest.approx(0.03)
    assert axial_condition(0.2, 14.0, "w", GEOM)[0] == pytest.approx(14.0 / GEOM.vtip)
    theta = collective_pitch(GEOM, ct_hover_target=0.02, theta0=None)
    assert theta > 0.0
    assert collective_pitch(GEOM, ct_hover_target=None, theta0=0.12) == pytest.approx(0.12)
    with pytest.raises(ValueError):
        collective_pitch(GEOM, ct_hover_target=0.02, theta0=0.12)


@pytest.mark.parametrize("model", ["uniform", "coleman_simple", "coleman_feingold", "drees"])
def test_coefficient_decompositions_all_models(model):
    """Verifies force and moment decompositions across all inflow models."""
    result = coefficients(0.2, 0.0, 0.12, GEOM, model)
    assert result["CQ"] == pytest.approx(result["CQi"] + result["CQ0"])
    assert result["CH"] == pytest.approx(result["CHi"] + result["CH0"])
    assert set(OUTPUTS).issubset(set(result.keys()))


def test_uniform_symmetry_and_coleman_gradient():
    """Verifies lateral symmetry for uniform inflow and non-zero moments for Coleman inflow."""
    uniform = coefficients(0.2, 0.0, 0.12, GEOM, "uniform")
    coleman = coefficients(0.2, 0.0, 0.12, GEOM, "coleman_feingold")
    assert uniform["CY"] == pytest.approx(0.0, abs=1e-14)
    assert uniform["CMy"] == pytest.approx(0.0, abs=1e-14)
    assert coleman["CY"] != pytest.approx(0.0, abs=1e-14)
    assert coleman["CMy"] > 0.0
    assert coleman["CMx"] < 0.0


def test_linear_twist_recovers_constant_when_delta_is_zero():
    """Ensures linear twist model reduces exactly to constant pitch when delta_twist is zero."""
    pitch_const = BladePitch("constant", 0.12, 0.12, 0.12, 0.20)
    pitch_linear = BladePitch("linear_twist", 0.12, 0.12, 0.12, 0.20)

    res_const = coefficients(0.3, 0.0, pitch_const, GEOM, "drees")
    res_linear = coefficients(0.3, 0.0, pitch_linear, GEOM, "drees")

    for k in OUTPUTS:
        assert res_linear[k] == pytest.approx(res_const[k], abs=1e-14)


def test_linear_twist_with_hover_target():
    """Validates that linear twist trims to the exact hover target CT."""
    target_ct = 0.022
    pitch = resolve_pitch(
        GEOM,
        pitch_mode="linear_twist",
        ct_hover_target=target_ct,
        theta_root_deg=14.0,
        theta_tip_deg=6.0,
    )
    res = coefficients(0.0, 0.0, pitch, GEOM, "uniform")
    assert res["CT"] == pytest.approx(target_ct, rel=1e-7)


def test_explicit_aerodynamic_models_are_the_defaults():
    """Verifies the default aerodynamic model selectors."""
    assert PROFILE_DRAG_MODEL == "numerical_vectorial"
    assert INDUCED_TORQUE_MODEL == "energy_balance"
    assert K_IND == pytest.approx(1.15)


def test_energy_balance_cqi_is_shaft_torque_not_power():
    """Verifies that energy-balance CQi excludes translational work mu*CHi."""
    mu, mu_z = 0.3, 0.03
    result = coefficients(mu, mu_z, 0.12, GEOM, "uniform")
    power_form = K_IND * result["lambda_i"] * result["CT"] + mu_z * result["CT"]
    assert result["CQi"] == pytest.approx(power_form - mu * result["CHi"])
    assert result["CQi"] < power_form


def test_aerodynamic_model_selectors_are_strict():
    """Rejects removed generic and superseded selector names instead of aliasing them."""
    for legacy in ("complete", "simple_bet", "analytical_bet", "numerical_profile"):
        with pytest.raises(ValueError, match="profile_drag_model"):
            profile_drag_coefficients(0.2, 0.0, GEOM, profile_drag_model=legacy)
    with pytest.raises(ValueError, match="induced_torque_model"):
        coefficients(
            0.2,
            0.0,
            0.12,
            GEOM,
            "uniform",
            induced_torque_model="complete",
        )
    with pytest.raises(ValueError, match="induced_torque_model"):
        coefficients(
            0.2,
            0.0,
            0.12,
            GEOM,
            "uniform",
            induced_torque_model="simple_bet",
        )


def test_three_profile_drag_paths_are_selectable():
    """Checks all three explicit profile-drag formulations and force decompositions."""
    tangential = coefficients(
        0.25, 0.01, 0.12, GEOM, "uniform",
        profile_drag_model="analytical_tangential",
        induced_torque_model="analytical_bet",
    )
    vectorial = coefficients(
        0.25, 0.01, 0.12, GEOM, "uniform",
        profile_drag_model="analytical_vectorial",
        induced_torque_model="analytical_bet",
    )
    numerical = coefficients(
        0.25, 0.01, 0.12, GEOM, "uniform",
        profile_drag_model="numerical_vectorial",
        induced_torque_model="energy_balance",
    )
    for result in (tangential, vectorial, numerical):
        assert result["CQ"] == pytest.approx(result["CQi"] + result["CQ0"])
        assert result["CH"] == pytest.approx(result["CHi"] + result["CH0"])
        assert "CT0" in result
    assert tangential["CT0"] == pytest.approx(0.0)
    assert vectorial["CT0"] < 0.0
    assert numerical["CT0"] < 0.0
    assert tangential["CH0"] < vectorial["CH0"]


def test_csv_is_named_zbet_and_separate_model_csvs(tmp_path):
    """Ensures CSV outputs include consolidated zBET.csv and separate zBET_<model>.csv files."""
    df = run_sweep(GEOM, np.array([0.0, 0.1]), "alpha", [0.0], 0.12, inflow_models=INFLOW_MODELS)
    res_path = save_results(df, tmp_path)
    assert res_path.name == "zBET.csv"
    assert (tmp_path / "zBET.csv").is_file()

    for model in INFLOW_MODELS:
        model_csv = tmp_path / f"zBET_{model}.csv"
        assert model_csv.is_file()
        m_df = pd.read_csv(model_csv)
        assert "model" in m_df.columns
        assert (m_df["model"] == model).all()
        assert "CPair" in m_df.columns
        assert "sigma_ref" in m_df.columns
        assert "sigma_geom" in m_df.columns
        assert "profile_drag_model" in m_df.columns
        assert "induced_torque_model" in m_df.columns
        assert (m_df["profile_drag_model"] == PROFILE_DRAG_MODEL).all()
        assert (m_df["induced_torque_model"] == INDUCED_TORQUE_MODEL).all()

    assert (tmp_path / "zBET_coleman.csv").is_file()


def test_solidity_input_equals_output_exactly():
    """Guarantees that input solidity (sigma_ref or sigma_geom) matches the corresponding output metric exactly."""
    val_ref = 0.20
    sol_ref = resolve_solidity("sigma_ref", sigma_ref=val_ref, root_cutout=0.20)
    assert sol_ref.sigma_ref == val_ref

    val_geom = 0.16
    sol_geom = resolve_solidity("sigma_geom", sigma_geom=val_geom, root_cutout=0.20)
    assert sol_geom.sigma_geom == val_geom


def test_sweep_and_separate_model_plots(tmp_path):
    """Validates advance ratio sweep execution and plot file generation."""
    df = run_sweep(GEOM, np.array([0.0, 0.2]), "alpha", [0.0, 10.0], 0.12, inflow_models=INFLOW_MODELS)
    assert len(df) == 16
    assert set(OUTPUTS).issubset(df.columns)
    assert set(df["model"]) == set(INFLOW_MODELS)
    assert "sigma_ref" in df.columns
    assert "sigma_geom" in df.columns
    assert "sigma_thrust" in df.columns
    assert "FoM" in df.columns
    assert "L_D_eff" in df.columns
    assert "CPair" in df.columns
    plot_results(df, tmp_path)
    pngs = list(tmp_path.glob("*.png"))
    # One plot per configured output and inflow model, plus Coleman aliases.
    assert len(pngs) >= len(OUTPUTS) * len(INFLOW_MODELS)
    assert any("cpair" in p.name for p in pngs)
    assert any("l_d_eff" in p.name for p in pngs)


def test_tip_loss_modes():
    """Validates blade tip loss modes (none, fixed, sissingh)."""
    # 1. Mode 'none' (B = 1.0):
    g_none = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20), tip_loss_mode="none")
    assert g_none.b_factor() == 1.0

    # 2. Mode 'fixed' (B = 0.96):
    g_fixed = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20), tip_loss_mode="fixed", tip_loss_b=0.96)
    assert g_fixed.b_factor() == pytest.approx(0.96)

    # 3. Mode 'sissingh' with CT = 0.022 and 4 blades: B = 1 - sqrt(2*0.022)/4
    g_siss = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20), tip_loss_mode="sissingh")
    b_siss = g_siss.b_factor(ct=0.022)
    expected_b = 1.0 - math.sqrt(2.0 * 0.022) / 4
    assert b_siss == pytest.approx(expected_b, rel=1e-5)

    # Physical effect: with tip loss enabled, lift-generated thrust decreases.
    res_none = coefficients(0.0, 0.0, 0.15, g_none, "uniform")
    res_fixed = coefficients(0.0, 0.0, 0.15, g_fixed, "uniform")
    assert res_fixed["CT"] < res_none["CT"]

    # Profile drag still acts over the physical blade to x=1 and must not be
    # truncated at the effective lift radius B.
    for profile_model in (
        "analytical_tangential",
        "analytical_vectorial",
        "numerical_vectorial",
    ):
        p_none = profile_drag_coefficients(0.3, 0.02, g_none, profile_model)
        p_fixed = profile_drag_coefficients(0.3, 0.02, g_fixed, profile_model)
        assert p_fixed == pytest.approx(p_none, rel=1e-12, abs=1e-14)


def test_prandtl_glauert_compressibility():
    """Validates Prandtl-Glauert compressibility correction."""
    # Disabled (incompressible baseline):
    g_incomp = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20), use_prandtl_glauert=False)
    assert g_incomp.lift_slope(mu=0.3) == pytest.approx(5.7)

    # Enabled:
    g_comp = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20), use_prandtl_glauert=True)
    a_comp = g_comp.lift_slope(mu=0.3)
    assert a_comp > 5.7
    res_incomp = coefficients(0.3, 0.0, 0.12, g_incomp, "uniform")
    res_comp = coefficients(0.3, 0.0, 0.12, g_comp, "uniform")
    assert res_comp["CT"] > res_incomp["CT"]


def test_hover_trim_collective_with_twist_preserves_delta_twist():
    """Ensures collective hover trim strictly preserves total twist delta_twist = const under linear twist."""
    geom = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20))
    th_root_deg = 15.0
    th_tip_deg = 3.0
    delta_twist_deg = th_tip_deg - th_root_deg  # -12.0 deg

    _, pitch_trimmed = trim_hover(
        geom,
        hover_trim_mode="collective",
        pitch_mode="linear_twist",
        ct_hover_target=0.022,
        theta_root_deg=th_root_deg,
        theta_tip_deg=th_tip_deg,
    )
    trimmed_delta_deg = math.degrees(pitch_trimmed.theta_tip - pitch_trimmed.theta_root)
    assert trimmed_delta_deg == pytest.approx(delta_twist_deg, abs=1e-12)

    delta_root = math.degrees(pitch_trimmed.theta_root) - th_root_deg
    delta_tip = math.degrees(pitch_trimmed.theta_tip) - th_tip_deg
    assert delta_root == pytest.approx(delta_tip, abs=1e-12)

    res_hover = coefficients(0.0, 0.0, pitch_trimmed, geom, "uniform")
    assert res_hover["CT"] == pytest.approx(0.022, rel=1e-4)


def test_hover_trim_thrust_in_newtons():
    """Validates hover thrust target input directly in Newtons."""
    geom = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20), rho=1.225)
    thrust_target_n = 3000.0

    _, pitch = trim_hover(
        geom,
        hover_trim_mode="collective",
        pitch_mode="constant",
        thrust_hover_n=thrust_target_n,
    )
    res = coefficients(0.0, 0.0, pitch, geom, "uniform")
    assert res["T_N"] == pytest.approx(thrust_target_n, rel=1e-4)


def test_hover_trim_rpm_mode():
    """Validates hover trim by RPM mode to achieve target thrust with fixed pitch."""
    geom = Geometry(1000.0, 1.4, 5.7, 0.20, 0.016, resolve_solidity("sigma_ref", sigma_ref=0.20), rho=1.225)
    thrust_target_n = 2500.0
    pitch_deg = 14.0

    geom_trimmed, pitch_fixed = trim_hover(
        geom,
        hover_trim_mode="rpm",
        pitch_mode="constant",
        thrust_hover_n=thrust_target_n,
        theta0=math.radians(pitch_deg),
    )
    assert math.degrees(pitch_fixed.theta0) == pytest.approx(pitch_deg)
    res = coefficients(0.0, 0.0, pitch_fixed, geom_trimmed, "uniform")
    assert res["T_N"] == pytest.approx(thrust_target_n, rel=1e-4)
