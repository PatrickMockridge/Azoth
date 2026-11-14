// The DEXPI 2.0 Process writer, for the `azoth dexpi export --profile process` oracle.
//
// **The writer is the oracle, and only the writer.** NeqSim has no native DEXPI 2.0 Process
// reader - `Dexpi20ProcessModelPackageReader` extracts and SHA-checks XML from a zip, it does
// not rebuild a `ProcessModel` - so an import cannot be checked against this probe. What this
// probe pins is what a `ProcessSystem` becomes when the reviewed writer serialises it: the
// steps, their types, the streams between them, and the quantities kg/h, bar absolute and °C.
//
// The topology is the one `crates/azoth-dexpi/tests/fixtures/oracle_simple.toml` mirrors:
// a feed, a separator, a compressor, a cooler, a throttling valve, and two product streams.
//
//     javac -proc:none -cp neqsim-f0c7436.jar Dexpi20ProcessProbe.java
//     java -cp .:neqsim-f0c7436.jar Dexpi20ProcessProbe > captures/dexpi_process_simple.tsv
//     java -cp .:neqsim-f0c7436.jar Dexpi20ProcessProbe simple \
//       --xml ../dexpi/neqsim/process_simple.dexpi.xml
//
// The default run writes no file and prints only the inventory, so `tools/oracle_sweep.py` is
// safe to run it like every other probe.

import java.io.File;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import neqsim.process.equipment.compressor.Compressor;
import neqsim.process.equipment.heatexchanger.Cooler;
import neqsim.process.equipment.separator.Separator;
import neqsim.process.equipment.stream.Stream;
import neqsim.process.equipment.valve.ThrottlingValve;
import neqsim.process.processmodel.ProcessSystem;
import neqsim.process.processmodel.dexpi.Dexpi20ProcessModelWriter;
import neqsim.process.processmodel.dexpi.Dexpi20SemanticValidator;
import neqsim.thermo.system.SystemSrkEos;

public final class Dexpi20ProcessProbe {
  private Dexpi20ProcessProbe() {
  }

  public static void main(String[] args) throws Exception {
    String caseName = args.length > 0 ? args[0] : "simple";
    String xmlPath = null;
    for (int index = 1; index < args.length; index++) {
      if ("--xml".equals(args[index]) && index + 1 < args.length) {
        xmlPath = args[++index];
      }
    }
    if (!"simple".equals(caseName)) {
      throw new IllegalArgumentException("unknown case: " + caseName);
    }

    ProcessSystem process = new ProcessSystem();
    process.setName("DEXPI process oracle");
    buildSimple(process);
    process.run();

    File file = xmlPath == null ? File.createTempFile("dexpi-process-", ".xml") : new File(xmlPath);
    Dexpi20ProcessModelWriter.write(process, file);
    System.out.println("# DEXPI 2.0 Process export inventory, pinned NeqSim");
    inventory(file.toPath());
  }

  /** The topology `oracle_simple.toml` mirrors: feed, separate, compress, cool, let down. */
  private static void buildSimple(ProcessSystem process) {
    SystemSrkEos fluid = new SystemSrkEos(303.15, 50.0);
    fluid.addComponent("methane", 0.8);
    fluid.addComponent("n-butane", 0.2);
    fluid.setMixingRule("classic");

    Stream feed = new Stream("FEED", fluid);
    feed.setFlowRate(1000.0, "kg/hr");
    feed.setTemperature(30.0, "C");
    feed.setPressure(50.0, "bara");

    Separator separator = new Separator("V-101", feed);
    Compressor compressor = new Compressor("K-101", separator.getGasOutStream());
    compressor.setOutletPressure(80.0);
    Cooler cooler = new Cooler("E-101", compressor.getOutletStream());
    cooler.setOutletTemperature(30.0, "C");
    ThrottlingValve valve = new ThrottlingValve("FV-101", cooler.getOutletStream());
    valve.setOutletPressure(10.0);
    Stream gas = new Stream("GAS", valve.getOutletStream());
    Stream oil = new Stream("OIL", separator.getLiquidOutStream());

    process.add(feed);
    process.add(separator);
    process.add(compressor);
    process.add(cooler);
    process.add(valve);
    process.add(gas);
    process.add(oil);
  }

  private static void inventory(Path path) throws Exception {
    Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(path.toFile());
    Element root = document.getDocumentElement();
    System.out.println("root.name=" + root.getAttribute("name"));
    for (Element element : elements(root, "Import")) {
      System.out.println("import." + element.getAttribute("prefix") + "=" + element.getAttribute("source"));
    }

    Map<String, String> portStep = new LinkedHashMap<>();
    Map<String, Element> steps = new LinkedHashMap<>();
    for (Element object : objects(root)) {
      String type = object.getAttribute("type");
      if (!type.startsWith("Process/Process.") || type.equals("Process/Process.MaterialPort")
          || type.equals("Process/Process.Stream")) {
        continue;
      }
      String identifier = data(object, "Identifier");
      if (identifier == null) {
        continue;
      }
      steps.put(identifier, object);
      System.out.println("step." + identifier + ".type=" + type);
      Element ports = component(object, "Ports");
      if (ports == null) {
        continue;
      }
      for (Element port : elements(ports, "Object")) {
        portStep.put(port.getAttribute("id"), identifier);
      }
    }

    int streams = 0;
    for (Element object : objects(root)) {
      if (!"Process/Process.Stream".equals(object.getAttribute("type"))) {
        continue;
      }
      streams++;
      String label = data(object, "Label");
      String source = firstReference(object, "Source");
      String target = firstReference(object, "Target");
      System.out.println("stream." + label + ".source="
          + portStep.getOrDefault(strip(source), "?"));
      System.out.println("stream." + label + ".target="
          + portStep.getOrDefault(strip(target), "?"));
      quantity(object, label, "MassFlow");
      quantity(object, label, "Pressure");
      quantity(object, label, "Temperature");
    }
    System.out.println("step.count=" + steps.size());
    System.out.println("stream.count=" + streams);

    Dexpi20SemanticValidator.ValidationReport report = Dexpi20SemanticValidator.validate(path);
    System.out.println("semantic.errors=" + report.getErrors().size());
    System.out.println("semantic.warnings=" + report.getWarnings().size());
  }

  private static void quantity(Element stream, String label, String property) {
    Element component = component(stream, property);
    if (component == null) {
      return;
    }
    for (Element qualified : elements(component, "Object")) {
      Element data = first(qualified, "Data");
      if (data == null) {
        continue;
      }
      Element aggregated = first(data, "AggregatedDataValue");
      if (aggregated == null) {
        continue;
      }
      for (Element entry : elements(aggregated, "Data")) {
        String name = entry.getAttribute("property");
        if ("Unit".equals(name)) {
          Element reference = first(entry, "DataReference");
          if (reference != null) {
            System.out.println("stream." + label + "." + property + ".unit=" + reference.getAttribute("data"));
          }
        } else if ("Value".equals(name)) {
          Element number = first(entry, "Double");
          if (number != null) {
            System.out.println("stream." + label + "." + property + "=" + number.getTextContent().trim());
          }
        }
      }
    }
  }

  private static Iterable<Element> objects(Element root) {
    NodeList all = root.getElementsByTagName("Object");
    java.util.List<Element> found = new java.util.ArrayList<>();
    for (int index = 0; index < all.getLength(); index++) {
      found.add((Element) all.item(index));
    }
    return found;
  }

  private static Iterable<Element> elements(Element parent, String tag) {
    NodeList all = parent.getElementsByTagName(tag);
    java.util.List<Element> found = new java.util.ArrayList<>();
    for (int index = 0; index < all.getLength(); index++) {
      Node node = all.item(index);
      if (node.getParentNode() == parent) {
        found.add((Element) node);
      }
    }
    return found;
  }

  private static Element first(Element parent, String tag) {
    for (Element element : elements(parent, tag)) {
      return element;
    }
    return null;
  }

  private static Element component(Element parent, String property) {
    for (Element element : elements(parent, "Components")) {
      if (property.equals(element.getAttribute("property"))) {
        return element;
      }
    }
    return null;
  }

  private static String data(Element parent, String property) {
    for (Element element : elements(parent, "Data")) {
      if (!property.equals(element.getAttribute("property"))) {
        continue;
      }
      Element value = first(element, "String");
      return value == null ? null : value.getTextContent().trim();
    }
    return null;
  }

  private static String firstReference(Element parent, String property) {
    for (Element element : elements(parent, "References")) {
      if (property.equals(element.getAttribute("property"))) {
        return element.getAttribute("objects").trim().split("\\s+")[0];
      }
    }
    return "";
  }

  private static String strip(String value) {
    return value.startsWith("#") ? value.substring(1) : value;
  }
}
