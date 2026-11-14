//! The writer is idempotent, and the reader recovers what the writer wrote.

mod common;

use azoth_dexpi::export_process;
use azoth_dexpi::process::export::ExportOptions;
use azoth_dexpi::reader::read;
use azoth_dexpi::writer::write;

#[test]
fn write_read_write_is_byte_identical() {
    let palette = common::palette();
    let flowsheet = common::demo();
    let options = ExportOptions {
        allow_lossy: true,
        project_recycles: true,
        azoth_extensions: true,
        ..ExportOptions::default()
    };
    let first = export_process(&flowsheet, &palette, &options)
        .expect("the demo exports")
        .xml;
    let document = read(&first).expect("the export reads").document;
    let second = write(&document).expect("the document writes");
    assert_eq!(first, second);
}

#[test]
fn the_reader_preserves_text_order_and_escapes() {
    let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<Model name="m" uri="urn:x">
    <Import prefix="Core" source="https://data.dexpi.org/models/2.0.0/Core.xml"/>
    <Object id="A" type="Core/EngineeringModel">
        <Data property="Label"><String>a &amp; b &lt; c</String></Data>
    </Object>
</Model>
"#;
    let document = read(xml).expect("the document reads").document;
    let object = document
        .model
        .children
        .iter()
        .find_map(|child| match child {
            azoth_dexpi::Child::Object(object) => Some(object),
            _ => None,
        })
        .expect("the object reads");
    assert_eq!(object.data_string("Label"), Some("a & b < c"));
    let rendered = write(&document).expect("the document writes");
    assert!(rendered.contains("a &amp; b &lt; c"));
}
