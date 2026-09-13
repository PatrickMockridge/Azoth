//! The shape of a calculation's specification, as the code sees it.
//!
//! `tools/gen_registry.py` turns the YAML under `specs/calcs/` into static tables
//! built from these types, and every calculation reads its own bounds and test
//! cases from those tables rather than restating them. That is what makes a spec
//! authoritative at runtime instead of merely descriptive.
//!
//! # Why these live in `azoth-core` and not in the generated file
//!
//! They were generated into `crates/azoth-hydraulics/src/spec_gen.rs`, which was
//! fine while there was exactly one namespace. With two, each namespace crate would
//! generate its *own* `CalcSpec`, `TestCase` and `SpecCheck` - distinct types with
//! the same names, one per crate. Anything taking `&CalcSpec` would then accept only
//! the copy from its own crate, so a shared helper or a cross-namespace list would
//! not compile, and the failure would look like a confusing mismatch of otherwise
//! identical structs.
//!
//! Declaring them once here means every namespace's generated tables are built from
//! the same types, and the only thing a namespace crate generates is its own data.
//! It also puts them where the rest of the shared vocabulary lives, which is what
//! this crate is for.
//!
//! # What is checked
//!
//! The generated tables are asserted against the specs by the contract tests, and
//! the `FIELDS` lists they carry are asserted against both the Python result
//! dataclasses and the Rust result structs. A spec edit that changes a bound
//! changes behaviour in both languages with no second edit; one that changes a
//! *name* fails a test rather than producing two implementations that agree
//! numerically and disagree about what to call the answer.

use crate::range::RangeCheck;

/// Solver configuration for an implicit calculation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverSpec {
    /// The scheme to solve with, as the spec spells it.
    ///
    /// A string rather than [`crate::solver::SolverKind`] because it is a direct
    /// transcription of the spec, and because the generated table is data rather
    /// than behaviour: the calc parses it with `SolverKind::parse`, which rejects
    /// an unknown name rather than defaulting to one that would run a scheme the
    /// spec did not ask for. `test_solver_contract.py` holds this vocabulary to the
    /// schema and to both implementations.
    pub kind: &'static str,
    /// Stopping tolerance.
    pub tolerance: f64,
    /// Iteration cap.
    pub max_iterations: u32,
    /// Starting value.
    ///
    /// `None` for a scheme that does not iterate from a declared start -
    /// `cubic_roots` forms the roots analytically and only *polishes* them, so a
    /// starting point is not part of its scheme. The schema requires this key
    /// conditionally for exactly that reason.
    pub initial_guess: Option<f64>,
    /// `absolute` or `relative`; see the spec schema for what each means.
    pub convergence: &'static str,
}

/// One test case from a spec's `tests` list, plus the worked example.
#[derive(Debug, Clone, Copy)]
pub struct TestCase {
    /// Test id, unique within the calc.
    pub id: &'static str,
    /// `worked_example`, `reference`, or `property`.
    pub kind: &'static str,
    /// Which invariant a `property` test checks.
    pub property: Option<&'static str>,
    /// `active` or `skipped`.
    pub status: &'static str,
    /// Why a skipped test does not run. Never `None` when status is `skipped`.
    pub skip_reason: Option<&'static str>,
    /// Relative tolerance for comparison.
    pub tolerance: f64,
    /// Scalar inputs, by name.
    pub numbers: &'static [(&'static str, f64)],
    /// List-valued inputs, by name - a list of fitting ids.
    pub lists: &'static [(&'static str, &'static [&'static str])],
    /// Vector-valued inputs, by name: one number per component, or a composition.
    ///
    /// Separate from [`Self::lists`] because the two are different things that both
    /// happen to be JSON arrays. A list of fitting ids is *identifiers*, which a
    /// model resolves against a data file; a vector is *numbers*, which it does
    /// arithmetic on. Collecting them with one rule gave a composition the type of a
    /// list of names.
    pub vectors: &'static [(&'static str, &'static [f64])],
    /// Matrix-valued inputs, by name, flattened row-major.
    ///
    /// The dimension is not carried: every matrix this registry declares is `N x N`
    /// for the same `N` its vectors have, so a consumer reshapes against the vector
    /// length rather than against a second copy of the same number.
    pub matrices: &'static [(&'static str, &'static [f64])],
    /// Expected scalar outputs, by name.
    pub expected: &'static [(&'static str, f64)],
    /// Expected vector outputs, by name - a phase composition, or one K-value per
    /// component. Separate from [`Self::expected`] for the same reason
    /// [`Self::vectors`] is separate from [`Self::lists`].
    pub expected_vectors: &'static [(&'static str, &'static [f64])],
}

impl TestCase {
    /// Fetch a scalar input. `None` if the spec does not supply it.
    #[must_use]
    pub fn input(&self, name: &str) -> Option<f64> {
        self.numbers
            .iter()
            .find(|(k, _)| *k == name)
            .map(|(_, v)| *v)
    }

    /// Fetch a list-valued input.
    #[must_use]
    pub fn list(&self, name: &str) -> Option<&'static [&'static str]> {
        self.lists.iter().find(|(k, _)| *k == name).map(|(_, v)| *v)
    }

    /// Fetch a vector-valued input.
    #[must_use]
    pub fn vector(&self, name: &str) -> Option<&'static [f64]> {
        self.vectors
            .iter()
            .find(|(k, _)| *k == name)
            .map(|(_, v)| *v)
    }

    /// Fetch a matrix-valued input, flattened row-major.
    #[must_use]
    pub fn matrix(&self, name: &str) -> Option<&'static [f64]> {
        self.matrices
            .iter()
            .find(|(k, _)| *k == name)
            .map(|(_, v)| *v)
    }

    /// Fetch an expected scalar output. `None` if the test does not assert it.
    #[must_use]
    pub fn expected_value(&self, name: &str) -> Option<f64> {
        self.expected
            .iter()
            .find(|(k, _)| *k == name)
            .map(|(_, v)| *v)
    }

    /// Fetch an expected vector output. `None` if the test does not assert it.
    #[must_use]
    pub fn expected_vector(&self, name: &str) -> Option<&'static [f64]> {
        self.expected_vectors
            .iter()
            .find(|(k, _)| *k == name)
            .map(|(_, v)| *v)
    }

    /// Whether this test should actually run.
    #[must_use]
    pub fn is_active(&self) -> bool {
        self.status == "active"
    }
}

/// A range check plus the phase in which it can be evaluated.
#[derive(Debug, Clone, Copy)]
pub struct SpecCheck {
    /// True when the bounded quantity is one of the calc's declared inputs, so
    /// the check can run before the calculation.
    ///
    /// False for outputs and derived quantities, which only exist afterwards -
    /// and, for an optional input, may not exist at all.
    pub on_input: bool,
    /// The bound itself.
    pub check: RangeCheck,
}

/// Everything the Rust side needs to know about one calculation.
#[derive(Debug, Clone, Copy)]
pub struct CalcSpec {
    /// Spec id, e.g. `hydraulics.darcy_weisbach`.
    pub id: &'static str,
    /// `verified`, `unverified` or `source_needed`.
    pub verification: &'static str,
    /// Bounds, in spec order.
    pub checks: &'static [SpecCheck],
    /// Present only for implicit calculations.
    pub solver: Option<SolverSpec>,
    /// The worked example, as a runnable test case.
    pub worked_example: TestCase,
    /// The rest of the `tests` list, excluding the worked example.
    pub tests: &'static [TestCase],
}

impl CalcSpec {
    /// Checks evaluable from the inputs alone.
    pub fn input_checks(&self) -> impl Iterator<Item = &RangeCheck> {
        self.checks.iter().filter(|c| c.on_input).map(|c| &c.check)
    }

    /// Checks that need the calculation to have run first.
    pub fn derived_checks(&self) -> impl Iterator<Item = &RangeCheck> {
        self.checks.iter().filter(|c| !c.on_input).map(|c| &c.check)
    }

    /// Every check, in spec order.
    pub fn all_checks(&self) -> impl Iterator<Item = &RangeCheck> {
        self.checks.iter().map(|c| &c.check)
    }

    /// The worked example plus every other test, in spec order.
    pub fn all_tests(&self) -> impl Iterator<Item = &TestCase> {
        std::iter::once(&self.worked_example).chain(self.tests.iter())
    }
}

/// A model's bracketing rule, as the spec declares it.
///
/// Part of the algorithm rather than an implementation detail. Two implementations
/// that bracket differently find different roots when there is more than one, and the
/// difference is invisible in a converged answer until it is a different answer.
#[derive(Debug, Clone, Copy)]
pub struct ModelBracket {
    /// How the bracket is found, as the spec spells it.
    pub scheme: &'static str,
    /// Where the search starts, in the reduced variable the scheme names - so the
    /// bracket does not depend on which unit a caller used.
    pub lower: f64,
    /// Where the search stops.
    pub upper: f64,
    /// How many points the search evaluates. Stated because it quantises the bracket:
    /// two implementations with different step counts take different numbers of
    /// iterations to the same answer.
    pub steps: u32,
}

/// A model's procedure, as the spec declares it.
///
/// This is the whole reason a model is not a calc. A calc's spec fixes an equation
/// and both implementations evaluate it; a model's fixes a *loop*, and a loop that
/// differs even slightly between two implementations diverges. `scheme` is held to
/// the implementations by a contract test, the way `SolverKind` is.
///
/// # Why this is recursive
///
/// [`Self::inner`] is another `ModelAlgorithm`. A nested procedure has the same
/// shape as an outer one - a scheme, a convergence rule, a tolerance, a cap - so
/// describing it with a second set of types would be a second definition of what an
/// algorithm is, free to drift from the first. The PT flash runs a Rachford-Rice
/// solve every outer iteration, and this is what lets the spec say so without
/// inventing a vocabulary for it.
#[derive(Debug, Clone, Copy)]
pub struct ModelAlgorithm {
    /// The procedure's name, e.g. `saturation_pressure_bisection`.
    pub scheme: &'static str,
    /// `absolute` or `relative`.
    pub convergence: &'static str,
    /// Stopping tolerance.
    pub tolerance: f64,
    /// Iteration cap.
    pub max_iterations: u32,
    /// How the starting interval is found, when the scheme brackets anything.
    ///
    /// Optional, because an inner scheme need not bracket: Rachford-Rice bisects a
    /// bracket the K-values define analytically, so there is no search for one and a
    /// mandatory field would have to be filled with a bracket nobody evaluates.
    pub bracket: Option<ModelBracket>,
    /// How the iteration starts, where the scheme needs a starting point.
    ///
    /// The analogue of a solver's `initial_guess`, and in the spec for the same
    /// reason: two implementations that start from different guesses take different
    /// paths to the same answer, and a convergence claim that holds from one
    /// starting point need not hold from another.
    pub initialisation: Option<&'static str>,
    /// A nested procedure the outer one runs each iteration.
    pub inner: Option<&'static ModelAlgorithm>,
}

/// One model, as the generated table carries it.
///
/// Deliberately reuses [`SpecCheck`] and [`TestCase`] from the calc side: a model's
/// bounds are the same kind of bound and its cases are the same kind of case, and a
/// second set of types for them would be a second definition of what a bound is.
#[derive(Debug, Clone, Copy)]
pub struct ModelSpec {
    /// Dotted identifier, e.g. `eos.pure_saturation`.
    pub id: &'static str,
    /// `verified`, `unverified` or `source_needed`, as for a calc.
    pub verification: &'static str,
    /// `procedure` or `direct` - what this spec fixes.
    ///
    /// A **procedure** is a loop, and [`Self::algorithm`] is present and is the whole
    /// contract. A **direct** model is a computation over vectors with no iteration,
    /// which the model tree holds because the calc registry's inputs are scalars and
    /// a composition vector has nowhere else to go; it has no algorithm, and the
    /// schema forbids one rather than allowing a vacuous one.
    pub kind: &'static str,
    /// The procedure both implementations must run, for a `procedure`.
    pub algorithm: Option<&'static ModelAlgorithm>,
    /// Bounds, in the same shape a calc's `valid_range` produces.
    pub checks: &'static [SpecCheck],
    /// Cases both implementations must reproduce.
    pub cases: &'static [TestCase],
}

impl ModelSpec {
    /// Checks evaluable from the inputs alone.
    pub fn input_checks(&self) -> impl Iterator<Item = &RangeCheck> {
        self.checks.iter().filter(|c| c.on_input).map(|c| &c.check)
    }

    /// Checks that need a computed value.
    pub fn derived_checks(&self) -> impl Iterator<Item = &RangeCheck> {
        self.checks.iter().filter(|c| !c.on_input).map(|c| &c.check)
    }
}
