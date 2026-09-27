//! `azoth crosscheck`, end to end.
//!
//! The sidecar's whole contract is the wire, so these drive the built binary rather than calling
//! the module: a request document on stdin, one verdict per line on stdout. What is asserted is
//! the *doctrine* as much as the arithmetic — that a refusal returns no number, that a divergence
//! is reported as a finding rather than a failure, and that a state azoth is not qualified on is
//! refused rather than answered with something plausible.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use serde_json::{Value, json};

/// One file the binary wrote. Removed when the test that made it ends.
struct Fixture(PathBuf);

impl Fixture {
    fn write(name: &str, body: &str) -> Self {
        let path = std::env::temp_dir().join(format!("azoth-crosscheck-{name}.toml"));
        std::fs::write(&path, body).expect("the fixture should be writable");
        Self(path)
    }

    fn argument(&self) -> String {
        self.0.display().to_string()
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// A request as a value, so a test varies one field rather than one substring.
///
/// The state is the one `pt_flash`'s documented example pins: Peng-Robinson, 330 K, 25 bara,
/// methane and n-butane at 0.6 and 0.4, whose vapour fraction is `0.8422055475803881`.
fn state() -> Value {
    json!({
        "model": "PR",
        "temperature": {"value": 330.0, "unit": "K"},
        "pressure": {"value": 25.0, "unit": "bara"},
        "flashType": "TP",
        "components": {"methane": 0.6, "n-butane": 0.4},
        "mixingRule": "classic",
    })
}

/// The wire form of a request: one line.
fn line(request: &Value) -> String {
    format!("{request}\n")
}

/// The binary, run the way a caller runs it: every line in at once, every answer out.
fn run(args: &[&str], input: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_azoth"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("the built binary should run");
    child
        .stdin
        .take()
        .expect("a stdin pipe")
        .write_all(input.as_bytes())
        .expect("the sidecar is reading");
    child
        .wait_with_output()
        .expect("the sidecar should exit when its stdin closes")
}

/// Every answer the binary produced for one run.
fn answers(args: &[&str], input: &str) -> Vec<Value> {
    let output = run(args, input);
    let text = String::from_utf8(output.stdout).expect("utf-8");
    text.lines()
        .map(|one| serde_json::from_str(one).unwrap_or_else(|error| panic!("{error}\n{one}")))
        .collect()
}

/// The single answer to one request.
fn one(args: &[&str], request: &Value) -> Value {
    let mut all = answers(args, &line(request));
    assert_eq!(all.len(), 1, "one request is one answer");
    all.pop().expect("one answer")
}

/// The answer to one request, with no card.
fn answer_of(request: &Value) -> Value {
    one(&["crosscheck"], request)
}

/// The state with a caller's own vapour fraction beside it.
fn reporting(vapour_fraction: Value) -> Value {
    let mut request = state();
    request["reported"] = json!({"vapourFraction": vapour_fraction});
    request
}

#[test]
fn agreement_with_the_reported_value_passes_the_gate() {
    let document = answer_of(&reporting(json!(0.842_205_547_580_388_1)));
    assert_eq!(document["status"], "success");
    assert_eq!(document["qualityGate"]["verdict"], "passed");
    assert_eq!(document["autoValidation"]["overall"], "PASS");
    assert_eq!(
        document["data"]["vapourFraction"]["magnitude_si"],
        0.842_205_547_580_388_1
    );
    assert_eq!(document["data"]["phase"], "two_phase");
}

#[test]
fn a_tolerance_the_caller_sets_is_the_one_judged_against() {
    // A difference of 5.5e-6: inside the default, outside a tolerance the caller tightens.
    let mut loose = reporting(json!(0.8422));
    loose["tolerance"] = json!(1e-3);
    assert_eq!(answer_of(&loose)["autoValidation"]["overall"], "PASS");

    let mut strict = reporting(json!(0.8422));
    strict["tolerance"] = json!(1e-9);
    assert_eq!(answer_of(&strict)["autoValidation"]["overall"], "WARNING");
}

#[test]
fn a_disagreement_is_a_finding_and_not_a_verdict_on_either_number() {
    let document = answer_of(&reporting(json!(0.9)));
    assert_eq!(document["status"], "success");
    assert_eq!(document["qualityGate"]["verdict"], "failed");
    assert_eq!(document["autoValidation"]["overall"], "WARNING");

    let check = &document["autoValidation"]["checks"][0];
    assert_eq!(check["rule"], "cross-implementation");
    assert_eq!(check["status"], "FAIL");
    // Both numbers and the tolerance travel, so the reader can judge rather than trust.
    let message = check["message"].as_str().expect("a sentence");
    assert!(message.contains("0.842205547580"), "{message}");
    assert!(message.contains("0.900000000000"), "{message}");
    assert!(message.contains("tolerance 1e-6"), "{message}");

    // And the answer says what a divergence is and is not evidence of.
    assert!(
        document["provenance"]["limitations"]
            .as_array()
            .expect("limitations")
            .iter()
            .any(|item| item
                .as_str()
                .is_some_and(|text| text.contains("transcription"))),
        "{document}"
    );
}

#[test]
fn a_model_azoth_cannot_express_is_refused_rather_than_answered_by_a_cubic() {
    let mut request = state();
    request["model"] = json!("CPA");
    let document = answer_of(&request);
    assert_eq!(document["status"], "blocked");
    // No number at all: that is the whole point of a refusal.
    assert!(document["data"].is_null(), "{document}");
    assert_eq!(document["autoValidation"]["checks"][0]["rule"], "model");
    let message = document["validation"]["message"]
        .as_str()
        .expect("the refusal says why");
    assert!(message.contains("CPA"), "{message}");
    assert!(
        message.contains("plausible-looking wrong number"),
        "{message}"
    );
}

#[test]
fn a_substance_the_databank_lacks_is_refused_by_name() {
    let mut request = state();
    request["components"] = json!({"unobtainium": 1.0});
    let document = answer_of(&request);
    assert_eq!(document["status"], "blocked");
    assert!(document["data"].is_null());
    let message = document["validation"]["message"]
        .as_str()
        .expect("a reason");
    assert!(message.contains("unobtainium"), "{message}");
    // The refusal names where a name comes from, so it is actionable.
    assert!(message.contains("COMP.csv"), "{message}");
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused_with_its_total() {
    let mut request = state();
    request["components"] = json!({"methane": 0.6, "n-butane": 0.3});
    let document = answer_of(&request);
    assert_eq!(document["status"], "blocked");
    let message = document["validation"]["message"]
        .as_str()
        .expect("a reason");
    assert!(message.contains("0.899"), "{message}");
    // The refusal says what the other side does, so the difference in policy is visible.
    assert!(message.contains("normalises"), "{message}");
    assert_eq!(document["warnings"][0]["code"], "COMPOSITION_NORMALIZED");
}

#[test]
fn a_negative_flash_root_is_not_compared_as_a_fraction() {
    // 25 C and 50 bara is single-phase vapour for this mixture. azoth reports the Rachford-Rice
    // root anyway, as the negative-flash value, and it is outside [0, 1] — comparing it against a
    // reported fraction would be a difference between two regimes, not between two numbers.
    let request = json!({
        "model": "SRK",
        "temperature": {"value": 25.0, "unit": "C"},
        "pressure": {"value": 50.0, "unit": "bara"},
        "flashType": "TP",
        "components": {"methane": 0.85, "ethane": 0.10, "propane": 0.05},
        "mixingRule": "classic",
        "reported": {"vapourFraction": 0.9},
    });
    let document = answer_of(&request);
    assert_eq!(document["data"]["phase"], "all_vapour");
    let message = document["autoValidation"]["checks"][0]["message"]
        .as_str()
        .expect("a sentence");
    assert!(message.contains("different regimes"), "{message}");
    assert!(message.contains("not comparable"), "{message}");
    // Flagged as a phase uncertainty rather than as a transcription difference.
    assert!(
        document["warnings"]
            .as_array()
            .expect("warnings")
            .iter()
            .any(|warning| warning["code"] == "TWO_PHASE_UNCERTAINTY"),
        "{document}"
    );
}

#[test]
fn an_answer_with_no_reported_value_is_not_a_failed_comparison() {
    let document = answer_of(&state());
    assert_eq!(document["status"], "success");
    assert_eq!(document["qualityGate"]["verdict"], "passed");
    let message = document["autoValidation"]["checks"][0]["message"]
        .as_str()
        .expect("a sentence");
    assert!(
        message.contains("nothing was checked against it"),
        "{message}"
    );
}

#[test]
fn a_card_overrides_a_substance_and_the_answer_moves() {
    let without = answer_of(&state());
    let card = Fixture::write(
        "override",
        "schema_version = 2\n[keyholder]\nname = \"test\"\n\n\
         [components.methane.Tc]\nvalue = 150.0\nunit = \"K\"\n",
    );
    let with = one(&["crosscheck", "--keycard", &card.argument()], &state());
    assert_ne!(
        without["data"]["vapourFraction"]["magnitude_si"],
        with["data"]["vapourFraction"]["magnitude_si"],
        "a card is the only thing in this CLI that can change an answer for a state"
    );
}

#[test]
fn a_card_is_refused_when_it_is_loaded_and_not_when_a_value_is_reached() {
    // An unknown section is data nothing would ever read, so the process does not start.
    let card = Fixture::write(
        "bad-section",
        "schema_version = 2\n[keyholder]\nname = \"test\"\n\n\
         [componants.methane.Tc]\nvalue = 150.0\nunit = \"K\"\n",
    );
    let output = run(
        &["crosscheck", "--keycard", &card.argument()],
        &line(&state()),
    );
    assert_eq!(
        output.status.code(),
        Some(2),
        "the card is refused as a whole"
    );
    assert!(output.stdout.is_empty(), "nothing is answered");
    let message = String::from_utf8(output.stderr).expect("utf-8");
    assert!(message.contains("componants"), "{message}");
    assert!(message.contains("expected one of"), "{message}");
}

#[test]
fn an_incomplete_substance_is_refused_when_it_is_resolved() {
    // Vocabulary is refused at load; completeness only bites when the substance is used, and then
    // it names the parameter that is missing rather than defaulting it.
    let card = Fixture::write(
        "incomplete",
        "schema_version = 2\n[keyholder]\nname = \"test\"\n\n\
         [components.unobtainium.Tc]\nvalue = 300.0\nunit = \"K\"\n\n\
         [components.unobtainium.omega]\nvalue = 0.1\nunit = \"dimensionless\"\n",
    );
    let mut request = state();
    request["components"] = json!({"unobtainium": 1.0});
    let document = one(&["crosscheck", "--keycard", &card.argument()], &request);
    assert_eq!(document["status"], "blocked");
    let message = document["validation"]["message"]
        .as_str()
        .expect("a reason");
    assert!(message.contains("Pc"), "{message}");
    assert!(message.contains("inventing data"), "{message}");
}

#[test]
fn a_line_that_is_not_a_request_is_answered_and_the_process_keeps_reading() {
    let input = format!("not json at all\n{}", line(&state()));
    let output = run(&["crosscheck"], &input);
    assert!(output.status.success(), "a bad line is not fatal");
    let all = answers(&["crosscheck"], &input);
    assert_eq!(all.len(), 2, "one answer per line, including the bad one");
    assert_eq!(all[0]["status"], "error");
    assert_eq!(all[1]["status"], "success");
}

#[test]
fn a_field_the_contract_may_add_is_ignored_rather_than_refused() {
    // Their contract says new optional fields may be added to inputs at any time, so a sidecar
    // that refused an unknown one would break on the release after this one.
    let mut request = state();
    request["somethingAddedLater"] = json!({"a": 1});
    assert_eq!(answer_of(&request)["status"], "success");
}

#[test]
fn a_flash_type_azoth_has_a_kernel_for_but_does_not_wire_here_is_refused() {
    let mut request = state();
    request["flashType"] = json!("PH");
    let document = answer_of(&request);
    assert_eq!(document["status"], "blocked");
    let message = document["validation"]["message"]
        .as_str()
        .expect("a reason");
    // The refusal points at the kernel that would answer it, so the gap is legible.
    assert!(message.contains("eos.ph_flash"), "{message}");
}

#[test]
fn the_envelope_carries_what_is_not_claimed_about_it() {
    let document = answer_of(&state());
    assert_eq!(document["apiVersion"], "1.0");
    assert_eq!(document["tool"], "azothCrossCheck");
    assert_eq!(
        document["provenance"]["limitations"]
            .as_array()
            .expect("limitations")
            .len(),
        4,
        "every answer says what the sidecar is not"
    );
    let assumptions = document["provenance"]["assumptions"]
        .as_array()
        .expect("assumptions");
    assert!(
        assumptions.iter().any(|item| item
            .as_str()
            .is_some_and(|text| text.contains("classical kij"))),
        "{document}"
    );
    // The default tolerance is a claim, so it travels with the answer rather than staying in the
    // source.
    assert!(
        document["provenance"]["assumptions"]
            .as_array()
            .expect("assumptions")
            .iter()
            .any(|item| item.as_str().is_some_and(|text| text.contains("keycard"))),
        "an answer with no card should say the databank answered"
    );
}

/// The repository root, for reaching the committed NeqSim material the demo is held to.
fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The committed `runFlash` requests the upstream demo is: agree, refuse a model, refuse a
/// substance.
fn demo_requests() -> String {
    let path = root().join("validation/neqsim/crosscheck/requests.ndjson");
    std::fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

#[test]
fn the_demo_requests_agree_and_refuse_as_the_demo_says() {
    // The file is the artefact a reader runs, so it is held to its own description here rather
    // than left to rot: one agreement and two refusals, in that order.
    let all = answers(&["crosscheck"], &demo_requests());
    assert_eq!(all.len(), 3, "one answer per request");
    assert_eq!(all[0]["autoValidation"]["overall"], "PASS", "{}", all[0]);
    assert_eq!(all[1]["status"], "blocked");
    assert_eq!(all[1]["autoValidation"]["checks"][0]["rule"], "model");
    assert_eq!(all[2]["status"], "blocked");
    assert_eq!(all[2]["autoValidation"]["checks"][0]["rule"], "components");
}

#[test]
fn the_demo_agrees_with_the_committed_neqsim_case() {
    // **This is what makes the agreement evidence rather than a round trip.** The reported value
    // below is not azoth's own: it is what NeqSim produced, committed beside the driver that
    // printed it and marked `verified`, and the demo asserts azoth lands on the same number
    // within that case's own tolerance. A cross-check that compared azoth against itself would
    // prove nothing and would pass whether or not the port was right.
    let path = root().join("validation/eos/methane_butane_flash_against_neqsim.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let case: Value = serde_json::from_str(&text).expect("a case document");

    // An unverified or source_needed case would make the demo circular, so the demo refuses to
    // ride one.
    assert_eq!(case["source"]["verification"], "verified", "{path:?}");
    assert!(
        case["source"]["attribution"]
            .as_str()
            .is_some_and(|text| text.contains("NeqSim")),
        "the case must name its upstream"
    );

    let neqsim = case["expected"]["vapour_fraction"]
        .as_f64()
        .expect("a vapour fraction");
    let tolerance = case["tolerance"].as_f64().expect("a tolerance");
    let inputs = &case["inputs"];

    let mut request = json!({
        "model": "PR",
        "temperature": {"value": inputs["T"], "unit": "K"},
        "pressure": {"value": inputs["P"], "unit": "Pa"},
        "flashType": "TP",
        "mixingRule": "classic",
        "reported": {"vapourFraction": neqsim},
    });
    let mut components = serde_json::Map::new();
    for (name, fraction) in inputs["components"]
        .as_array()
        .expect("names")
        .iter()
        .zip(inputs["z"].as_array().expect("fractions"))
    {
        components.insert(name.as_str().expect("a name").to_string(), fraction.clone());
    }
    request["components"] = Value::Object(components);

    let document = answer_of(&request);
    assert_eq!(document["status"], "success", "{document}");
    assert_eq!(
        document["qualityGate"]["verdict"], "passed",
        "azoth and NeqSim agree on this state: {document}"
    );

    // And the agreement is the case's, not the default tolerance's: measured at 3.3e-12.
    let ours = document["data"]["vapourFraction"]["magnitude_si"]
        .as_f64()
        .expect("a vapour fraction");
    let difference = (ours - neqsim).abs();
    assert!(
        difference <= tolerance,
        "azoth {ours} against NeqSim {neqsim} is {difference:e}, outside the case's {tolerance:e}"
    );
}
