//! Builders shared by the unit-operation tests.
//!
//! The equation-of-state crate has the same two builders in `crates/azoth-eos/tests/` and
//! they cannot be shared with these: `azoth-test-support` depends on `azoth-core` alone,
//! so it cannot construct a `Mixture` or an `IdealGasModel`. Putting them here rather
//! than in a `tests/common/` under every crate is the next best thing - **this is the
//! only copy outside `azoth-eos`**, and a third would be the point at which they belong
//! in a shared test crate that depends on the domains.

use azoth_core::spec::TestCase;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::mixture::{Component, Mixture};
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;

/// The mixture a spec case describes, from its four per-component vectors.
///
/// Panics rather than returning a `Result`: a case missing `Tc` is a defect in the spec,
/// and a test that silently skipped it would be a test that reports success for having
/// run nothing.
pub fn mixture_from_case(case: &TestCase) -> Mixture {
    let tc = case.vector("Tc").expect("the case declares Tc");
    let pc = case.vector("Pc").expect("the case declares Pc");
    let omega = case.vector("omega").expect("the case declares omega");
    let kij = case.matrix("kij").expect("the case declares kij");

    let components = (0..tc.len())
        .map(|i| {
            Component::new(kelvins(tc[i]), pascals(pc[i]), omega[i]).expect("a valid component")
        })
        .collect();
    Mixture::new(components, kij.to_vec()).expect("a valid mixture")
}

/// The ideal-gas datum a spec case describes.
///
/// Every unit operation except `process.splitter` declares one, because every one except
/// a splitter does an energy balance. A case for a model that takes none will not have
/// these inputs, and the caller is the one that knows which it is asking for.
pub fn ideal_gas_from_case(case: &TestCase) -> IdealGasModel {
    fn vector(case: &TestCase, name: &str) -> Vec<f64> {
        case.vector(name)
            .unwrap_or_else(|| panic!("case `{}` declares no {name}", case.id))
            .to_vec()
    }
    IdealGasModel {
        cp_a: vector(case, "cp_a"),
        cp_b: vector(case, "cp_b"),
        cp_c: vector(case, "cp_c"),
        cp_d: vector(case, "cp_d"),
        h_ref: vector(case, "h_ref"),
        s_ref: vector(case, "s_ref"),
        t_ref: kelvins(azoth_test_support::input(case, "T_ref")),
        p_ref: pascals(azoth_test_support::input(case, "P_ref")),
    }
}
