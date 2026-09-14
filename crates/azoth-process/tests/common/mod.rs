//! Builders shared by the unit-operation tests.
//!
//! The equation-of-state crate has the same two builders in `crates/azoth-eos/tests/` and
//! they cannot be shared with these: `azoth-test-support` depends on `azoth-core` alone,
//! so it cannot construct a `Mixture` or an `IdealGasModel`. Putting them here rather
//! than in a `tests/common/` under every crate is the next best thing - **this is the
//! only copy outside `azoth-eos`**, and a third would be the point at which they belong
//! in a shared test crate that depends on the domains.

use azoth_core::spec::TestCase;
use azoth_eos::databank;
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;

/// The mixture a spec case describes, resolved from its `components` list.
///
/// Panics rather than returning a `Result`: a case whose components do not resolve is a
/// defect in the spec, and a test that silently skipped it would be a test that reports
/// success for having run nothing.
pub fn mixture_from_case(case: &TestCase) -> Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names)
        .expect("the case's components resolve")
        .0
}

/// The ideal-gas datum a spec case describes.
///
/// Every unit operation except `process.splitter` declares one, because every one except
/// a splitter does an energy balance. A case for a model that takes none will not have
/// these inputs, and the caller is the one that knows which it is asking for.
pub fn ideal_gas_from_case(case: &TestCase) -> IdealGasModel {
    let names = case
        .list("components")
        .expect("the case declares components");
    databank::mixture_of(names)
        .expect("the case's components resolve")
        .1
}
