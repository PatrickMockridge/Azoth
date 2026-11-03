//! `process.separator_performance`: the curve refusals and the branch a case cannot reach.

use azoth_process::kernels::separator_performance::{
    CurveSettings, DistributionKind, DropletDistribution, Internals, LookupPoint, Orientation,
    PerformanceSettings, separator_performance,
};

fn distribution(cut: f64) -> DropletDistribution {
    DropletDistribution {
        kind: DistributionKind::RosinRammler,
        characteristic_diameter: cut,
        spread: 2.6,
    }
}

fn base() -> PerformanceSettings {
    PerformanceSettings {
        gas_density: 50.0,
        oil_density: 800.0,
        water_density: 1020.0,
        gas_viscosity: 1.5e-5,
        oil_viscosity: 2.0e-3,
        water_viscosity: 8.0e-4,
        gas_velocity: 1.0,
        vessel_diameter: 2.0,
        vessel_length: 8.0,
        orientation: Orientation::Horizontal,
        liquid_level_fraction: 0.5,
        gas_liquid: Some(distribution(100.0e-6)),
        gas_bubble: None,
        oil_in_water: None,
        water_in_oil: None,
        mist_eliminator: None,
        coalescer: None,
        include_gravity_section: true,
        liquid_in_gas_factor: 1.0,
        gas_carry_under_factor: 1.0,
        liquid_liquid_factor: 1.0,
        enhanced: false,
        inlet_pipe_diameter: 0.2,
        surface_tension: 0.025,
        oil_volume_fraction: 0.5,
        apply_turbulence_correction: true,
        inlet_device: azoth_process::separator_design::inlet_device::InletDevice::HalfPipe,
        design_k_factor: None,
        liquid_liquid_residence_time: None,
        mist_position: 0.90,
        weir_position: 0.60,
        override_bulk_efficiency: 0.0,
        override_dsd_multiplier: 0.0,
        oil_water_interfacial_tension: 0.030,
    }
}

#[test]
fn a_custom_curve_needs_its_table() {
    let mut settings = base();
    settings.mist_eliminator = Some(CurveSettings {
        internals: Internals::Custom,
        cut_diameter: None,
        sharpness: None,
        max_efficiency: None,
        points: None,
    });
    assert!(separator_performance(&settings).is_err());
}

#[test]
fn a_plate_pack_needs_its_cut_and_efficiency() {
    let mut settings = base();
    settings.mist_eliminator = Some(CurveSettings {
        internals: Internals::PlatePack,
        cut_diameter: None,
        sharpness: None,
        max_efficiency: Some(0.98),
        points: None,
    });
    assert!(separator_performance(&settings).is_err());
}

#[test]
fn a_lookup_table_on_a_wire_mesh_is_refused() {
    let mut settings = base();
    settings.mist_eliminator = Some(CurveSettings {
        internals: Internals::WireMesh,
        cut_diameter: None,
        sharpness: None,
        max_efficiency: None,
        points: Some(vec![
            LookupPoint {
                diameter: 1.0e-6,
                efficiency: 0.1,
            },
            LookupPoint {
                diameter: 1.0e-5,
                efficiency: 0.9,
            },
        ]),
    });
    assert!(separator_performance(&settings).is_err());
}

#[test]
fn no_gas_liquid_distribution_leaves_its_branch_unrun() {
    let mut settings = base();
    settings.gas_liquid = None;
    settings.water_in_oil = Some(distribution(100.0e-6));
    settings.oil_in_water = Some(distribution(80.0e-6));
    let outcome = separator_performance(&settings).expect("the liquid-liquid stages answer");

    assert_eq!(outcome.oil_in_gas, 0.0);
    assert_eq!(outcome.water_in_gas, 0.0);
    assert_eq!(outcome.gravity_cut_diameter, 0.0);
    assert_eq!(outcome.overall_gas_liquid_efficiency, 0.0);
    // The liquid-liquid stages still ran, which is the class's own bare-constructor state.
    assert!(outcome.water_in_oil > 0.0);
    assert!(outcome.oil_in_water > 0.0);
    assert!(outcome.liquid_liquid_gravity_efficiency > 0.0);
}

#[test]
fn the_three_phase_entrainment_is_written_to_both_liquids() {
    let mut settings = base();
    settings.mist_eliminator = Some(CurveSettings {
        internals: Internals::WireMesh,
        cut_diameter: None,
        sharpness: None,
        max_efficiency: None,
        points: None,
    });
    let outcome = separator_performance(&settings).expect("the default curve");
    assert_eq!(outcome.oil_in_gas, outcome.water_in_gas);
    assert!(outcome.oil_in_gas > 0.0);
}
