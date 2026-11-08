//! `BellDelawareMethod` - the shell-side helpers a shell-and-tube rating reads.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/BellDelawareMethod.java` at
//! NeqSim `f0c7436c69`. Every coefficient is the class's own literal; where a published
//! correlation would differ, the class's inline constant is what is carried.

/// The ideal crossflow coefficient for a tube bank, Zhukauskas' two layouts.
///
/// `triangular` selects the staggered branch; the square branch differs only in its
/// high-Reynolds coefficient (`0.033` against `0.021`) and its exponent (`0.84` against
/// `0.8`).
pub fn ideal_crossflow_htc(
    mass_flux: f64,
    tube_od: f64,
    mu: f64,
    cp: f64,
    k: f64,
    mu_wall: f64,
    triangular: bool,
) -> f64 {
    if mu <= 0.0 || k <= 0.0 || cp <= 0.0 || tube_od <= 0.0 {
        return 0.0;
    }
    let re = mass_flux * tube_od / mu;
    let pr = cp * mu / k;
    let pr_wall = cp * mu_wall / k;
    let viscosity_correction = if pr_wall > 0.0 {
        (pr / pr_wall).powf(0.25)
    } else {
        1.0
    };
    let nu = if triangular {
        if re < 100.0 {
            0.9 * re.powf(0.4) * pr.powf(0.36) * viscosity_correction
        } else if re < 1000.0 {
            0.52 * re.powf(0.5) * pr.powf(0.36) * viscosity_correction
        } else if re < 200_000.0 {
            0.27 * re.powf(0.63) * pr.powf(0.36) * viscosity_correction
        } else {
            0.033 * re.powf(0.8) * pr.powf(0.36) * viscosity_correction
        }
    } else if re < 100.0 {
        0.9 * re.powf(0.4) * pr.powf(0.36) * viscosity_correction
    } else if re < 1000.0 {
        0.52 * re.powf(0.5) * pr.powf(0.36) * viscosity_correction
    } else if re < 200_000.0 {
        0.27 * re.powf(0.63) * pr.powf(0.36) * viscosity_correction
    } else {
        0.021 * re.powf(0.84) * pr.powf(0.36) * viscosity_correction
    };
    nu * k / tube_od
}

/// Kern's shell-side coefficient: `0.36 (k/De) Re**0.55 Pr**(1/3) (mu/mu_w)**0.14`.
pub fn kern_shell_side_htc(
    mass_flux: f64,
    equivalent_diameter: f64,
    mu: f64,
    cp: f64,
    k: f64,
    mu_wall: f64,
) -> f64 {
    if mu <= 0.0 || k <= 0.0 || cp <= 0.0 || equivalent_diameter <= 0.0 {
        return 0.0;
    }
    let re = mass_flux * equivalent_diameter / mu;
    let pr = cp * mu / k;
    let viscosity_correction = if mu_wall > 0.0 {
        (mu / mu_wall).powf(0.14)
    } else {
        1.0
    };
    let nu = 0.36 * re.powf(0.55) * pr.powf(1.0 / 3.0) * viscosity_correction;
    nu * k / equivalent_diameter
}

/// Kern's shell-side pressure drop.
pub fn kern_shell_side_pressure_drop(
    mass_flux: f64,
    equivalent_diameter: f64,
    shell_diameter: f64,
    baffle_count: usize,
    rho: f64,
    mu: f64,
    mu_wall: f64,
) -> f64 {
    if mu <= 0.0 || rho <= 0.0 || equivalent_diameter <= 0.0 {
        return 0.0;
    }
    let re = mass_flux * equivalent_diameter / mu;
    let viscosity_correction = if mu_wall > 0.0 {
        (mu / mu_wall).powf(0.14)
    } else {
        1.0
    };
    let friction = if re < 10.0 {
        64.0 / re.max(0.1)
    } else if re < 1000.0 {
        5.6 * re.powf(-0.344)
    } else {
        0.228 * re.powf(-0.1574)
    };
    let velocity = mass_flux / rho;
    friction * shell_diameter * (baffle_count as f64 + 1.0) * rho * velocity * velocity
        / (2.0 * equivalent_diameter * viscosity_correction)
}

/// The shell-side equivalent diameter for the two layouts.
pub fn shell_equivalent_diameter(tube_od: f64, tube_pitch: f64, triangular: bool) -> f64 {
    if tube_od <= 0.0 || tube_pitch <= 0.0 {
        return 0.0;
    }
    if triangular {
        let flow_area = (3.0_f64).sqrt() / 4.0 * tube_pitch * tube_pitch
            - std::f64::consts::PI / 8.0 * tube_od * tube_od;
        let wetted_perimeter = std::f64::consts::PI / 2.0 * tube_od;
        4.0 * flow_area / wetted_perimeter
    } else {
        let flow_area = tube_pitch * tube_pitch - std::f64::consts::PI / 4.0 * tube_od * tube_od;
        let wetted_perimeter = std::f64::consts::PI * tube_od;
        4.0 * flow_area / wetted_perimeter
    }
}

/// The shell-side crossflow area at the bundle centreline.
pub fn crossflow_area(shell_id: f64, baffle_spacing: f64, tube_od: f64, tube_pitch: f64) -> f64 {
    if shell_id <= 0.0 || baffle_spacing <= 0.0 || tube_pitch <= 0.0 {
        return 0.0;
    }
    baffle_spacing * shell_id * (1.0 - tube_od / tube_pitch)
}

/// `J_c`, the baffle-cut correction factor: `0.55 + 0.72 F_c`.
pub fn jc(baffle_cut: f64) -> f64 {
    0.55 + 0.72 * crossflow_fraction(baffle_cut)
}

/// `J_l`, the baffle-leakage correction factor.
pub fn jl(
    tube_to_baffle_clearance: f64,
    shell_to_baffle_clearance: f64,
    crossflow_area: f64,
    tube_count: usize,
    tube_od: f64,
    baffle_spacing: f64,
) -> f64 {
    if crossflow_area <= 0.0 {
        return 0.7;
    }
    let tube_baffle_area =
        tube_count as f64 * std::f64::consts::PI * tube_od * tube_to_baffle_clearance;
    let shell_baffle_area = std::f64::consts::PI * baffle_spacing * shell_to_baffle_clearance;
    let total_leakage_ratio = (tube_baffle_area + shell_baffle_area) / crossflow_area;
    let shell_leakage_fraction =
        shell_baffle_area / (tube_baffle_area + shell_baffle_area).max(1e-10);
    let value = 0.44 * (1.0 - shell_leakage_fraction)
        + (1.0 - 0.44 * (1.0 - shell_leakage_fraction)) * (-2.2 * total_leakage_ratio).exp();
    value.clamp(0.2, 1.0)
}

/// `J_b`, the bundle-bypass correction factor.
pub fn jb(
    bypass_area: f64,
    crossflow_area: f64,
    has_sealing: bool,
    sealing_pairs: usize,
    tube_rows_crossflow: usize,
) -> f64 {
    if crossflow_area <= 0.0 {
        return 0.8;
    }
    let fbp = bypass_area / crossflow_area;
    if has_sealing && tube_rows_crossflow > 0 {
        let ratio = sealing_pairs as f64 / tube_rows_crossflow as f64;
        if ratio >= 0.5 {
            return 1.0;
        }
    }
    (-1.35 * fbp).exp().clamp(0.5, 1.0)
}

/// `J_s`, the unequal-spacing correction factor.
pub fn js(
    central_spacing: f64,
    inlet_spacing: f64,
    outlet_spacing: f64,
    baffle_count: usize,
) -> f64 {
    if central_spacing <= 0.0 || baffle_count < 1 {
        return 1.0;
    }
    let nb = baffle_count as f64;
    let inlet_ratio = central_spacing / inlet_spacing;
    let outlet_ratio = central_spacing / outlet_spacing;
    let numerator = nb - 1.0 + inlet_ratio.powf(0.6) + outlet_ratio.powf(0.6);
    let denominator = nb - 1.0 + inlet_ratio + outlet_ratio;
    if denominator <= 0.0 {
        return 1.0;
    }
    (numerator / denominator).clamp(0.5, 1.0)
}

/// `J_r`, the adverse-temperature-gradient correction factor.
pub fn jr(re: f64, tube_rows_crossflow: usize) -> f64 {
    if re >= 100.0 || tube_rows_crossflow <= 10 {
        return 1.0;
    }
    let laminar = (10.0 / tube_rows_crossflow as f64).powf(0.18);
    if re < 20.0 {
        return laminar.clamp(0.4, 1.0);
    }
    let fraction = (re - 20.0) / 80.0;
    let value = laminar + fraction * (1.0 - laminar);
    value.clamp(0.4, 1.0)
}

/// The fully corrected Bell-Delaware coefficient.
pub fn corrected_htc(ideal: f64, jc: f64, jl: f64, jb: f64, js: f64, jr: f64) -> f64 {
    ideal * jc * jl * jb * js * jr
}

/// `R_l`, the leakage correction for pressure drop.
pub fn rl(
    tube_to_baffle_clearance: f64,
    shell_to_baffle_clearance: f64,
    crossflow_area: f64,
    tube_count: usize,
    tube_od: f64,
    baffle_spacing: f64,
) -> f64 {
    if crossflow_area <= 0.0 {
        return 0.5;
    }
    let tube_baffle_area =
        tube_count as f64 * std::f64::consts::PI * tube_od * tube_to_baffle_clearance;
    let shell_baffle_area = std::f64::consts::PI * baffle_spacing * shell_to_baffle_clearance;
    let total_leakage_ratio = (tube_baffle_area + shell_baffle_area) / crossflow_area;
    let shell_leakage_fraction =
        shell_baffle_area / (tube_baffle_area + shell_baffle_area).max(1e-10);
    let value = (-1.33 * (1.0 + shell_leakage_fraction) * total_leakage_ratio.powf(0.8)).exp();
    value.clamp(0.1, 1.0)
}

/// `R_b`, the bypass correction for pressure drop.
pub fn rb(
    bypass_area: f64,
    crossflow_area: f64,
    has_sealing: bool,
    sealing_pairs: usize,
    tube_rows_crossflow: usize,
) -> f64 {
    if crossflow_area <= 0.0 {
        return 0.7;
    }
    let fbp = bypass_area / crossflow_area;
    if has_sealing && tube_rows_crossflow > 0 {
        let ratio = sealing_pairs as f64 / tube_rows_crossflow as f64;
        if ratio >= 0.5 {
            return 1.0;
        }
    }
    (-3.7 * fbp).exp().clamp(0.3, 1.0)
}

/// The ideal crossflow pressure drop for one baffle compartment.
pub fn ideal_crossflow_pressure_drop(
    tube_rows_crossflow: usize,
    mass_flux: f64,
    rho: f64,
    re: f64,
    triangular: bool,
) -> f64 {
    if rho <= 0.0 || mass_flux <= 0.0 {
        return 0.0;
    }
    let friction = if triangular {
        if re < 1000.0 {
            48.0 / re.max(0.1)
        } else {
            0.35 * re.powf(-0.15)
        }
    } else if re < 1000.0 {
        35.0 / re.max(0.1)
    } else {
        0.18 * re.powf(-0.15)
    };
    friction * tube_rows_crossflow as f64 * mass_flux * mass_flux / (2.0 * rho)
}

/// The window-zone pressure drop per baffle.
pub fn window_pressure_drop(
    window_area: f64,
    mass_flow: f64,
    rho: f64,
    tube_rows_window: usize,
) -> f64 {
    if window_area <= 0.0 || rho <= 0.0 {
        return 0.0;
    }
    let velocity = mass_flow / (rho * window_area);
    (2.0 + 0.6 * tube_rows_window as f64) * rho * velocity * velocity / 2.0
}

/// The fraction of tubes in the crossflow zone, linearly interpolated between the
/// class's two bounds.
pub fn crossflow_fraction(baffle_cut: f64) -> f64 {
    if baffle_cut <= 0.15 {
        return 0.9;
    }
    if baffle_cut >= 0.45 {
        return 0.3;
    }
    0.9 - (baffle_cut - 0.15) / (0.45 - 0.15) * 0.6
}

/// The number of tube rows in the crossflow zone.
pub fn tube_rows_crossflow(
    shell_id: f64,
    baffle_cut: f64,
    tube_pitch: f64,
    triangular: bool,
) -> usize {
    if shell_id <= 0.0 || tube_pitch <= 0.0 {
        return 1;
    }
    let distance_between_tips = shell_id * (1.0 - 2.0 * baffle_cut);
    let row_pitch = if triangular {
        tube_pitch * (3.0_f64).sqrt() / 2.0
    } else {
        tube_pitch
    };
    ((distance_between_tips / row_pitch) as usize).max(1)
}

/// The bypass area between the outermost tubes and the shell wall.
pub fn bypass_area(shell_id: f64, bundle_diameter: f64, baffle_spacing: f64) -> f64 {
    if shell_id <= 0.0 || bundle_diameter <= 0.0 || baffle_spacing <= 0.0 {
        return 0.0;
    }
    (shell_id - bundle_diameter) * baffle_spacing
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn js_is_one_when_the_three_spacings_agree() {
        assert!((js(0.2, 0.2, 0.2, 10) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn the_crossflow_fraction_clamps_at_both_bounds() {
        assert!((crossflow_fraction(0.1) - 0.9).abs() < 1e-15);
        assert!((crossflow_fraction(0.5) - 0.3).abs() < 1e-15);
        assert!((crossflow_fraction(0.3) - 0.6).abs() < 1e-15);
    }

    #[test]
    fn sealing_above_half_the_rows_returns_one() {
        assert!((jb(0.01, 0.02, true, 6, 10) - 1.0).abs() < 1e-15);
        assert!((rb(0.01, 0.02, true, 6, 10) - 1.0).abs() < 1e-15);
        assert!(jb(0.01, 0.02, false, 6, 10) < 1.0);
    }

    #[test]
    fn the_equivalent_diameter_is_positive_for_both_layouts() {
        assert!(shell_equivalent_diameter(0.025, 0.032, true) > 0.0);
        assert!(shell_equivalent_diameter(0.025, 0.032, false) > 0.0);
    }
}
