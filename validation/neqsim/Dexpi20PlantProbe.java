// The DEXPI 2.0 Plant writer, for the `azoth dexpi export --profile plant` oracle.
//
// The same small `ProcessSystem` the Process probe builds, serialised as a P&ID: tagged items
// with nozzles, and one piping network system whose segments join them. The probe prints the
// inventory the Rust test compares - item tag and type, nozzle sub-tags, and the internal
// edges between tagged items - from the document it wrote, not from the objects it built.
//
//     javac -proc:none -cp neqsim-f0c7436.jar Dexpi20PlantProbe.java
//     java -cp .:neqsim-f0c7436.jar Dexpi20PlantProbe > captures/dexpi_plant_simple.tsv
//     java -cp .:neqsim-f0c7436.jar Dexpi20PlantProbe --xml ../dexpi/neqsim/plant_simple.dexpi.xml
//
// NeqSim's writer draws boundary segments only when `Dexpi20PlantExportOptions` asks for
// off-page connectors; the default run does not, so the reference has the internal edges and
// this port's export has those plus its boundary connectors. The comparison is on the items
// and the internal edges, which is what both writers state the same way.

import java.io.File;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import neqsim.process.equipment.compressor.Compressor;
import neqsim.process.equipment.heatexchanger.Cooler;
import neqsim.process.equipment.separator.Separator;
import neqsim.process.equipment.stream.Stream;
import neqsim.process.processmodel.ProcessSystem;
import neqsim.process.processmodel.dexpi.Dexpi20SemanticValidator;
import neqsim.process.processmodel.dexpi.Dexpi20XmlWriter;
import neqsim.thermo.system.SystemSrkEos;

public final class Dexpi20PlantProbe {
  private Dexpi20PlantProbe() {
  }

  public static void main(String[] args) throws Exception {
    String xmlPath = null;
    for (int index = 0; index < args.length; index++) {
      if ("--xml".equals(args[index]) && index + 1 < args.length) {
        xmlPath = args[++index];
      }
    }

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
    Stream gas = new Stream("GAS", cooler.getOutletStream());
    Stream oil = new Stream("OIL", separator.getLiquidOutStream());

    ProcessSystem process = new ProcessSystem();
    process.setName("DEXPI plant oracle");
    process.add(feed);
    process.add(separator);
    process.add(compressor);
    process.add(cooler);
    process.add(gas);
    process.add(oil);

    File file = xmlPath == null ? File.createTempFile("dexpi-plant-", ".xml") : new File(xmlPath);
    Dexpi20XmlWriter.write(process, file);
    System.out.println("# DEXPI 2.0 Plant export inventory, pinned NeqSim");
    inventory(file.toPath());
  }

  private static void inventory(Path path) throws Exception {
    Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(path.toFile());
    Element root = document.getDocumentElement();
    System.out.println("root.name=" + root.getAttribute("name"));
    for (Element element : elements(root, "Import")) {
      System.out.println("import." + element.getAttribute("prefix") + "=" + element.getAttribute("source"));
    }

    Map<String, String> nodeItem = new LinkedHashMap<>();
    Map<String, String> itemTag = new LinkedHashMap<>();
    int items = 0;
    int nozzles = 0;
    for (Element object : objects(root)) {
      String type = object.getAttribute("type");
      if (!type.startsWith("Plant/ProcessEquipment.") || type.equals("Plant/ProcessEquipment.Nozzle")) {
        continue;
      }
      String tag = data(object, "TagName");
      if (tag == null) {
        continue;
      }
      items++;
      String id = object.getAttribute("id");
      itemTag.put(id, tag);
      System.out.println("item." + tag + ".type=" + type);
      Element nozzleList = component(object, "Nozzles");
      if (nozzleList == null) {
        continue;
      }
      for (Element nozzle : elementsOf(nozzleList, "Object")) {
        nozzles++;
        String subTag = data(nozzle, "SubTagName");
        System.out.println("item." + tag + ".nozzle=" + subTag);
        Element nodes = component(nozzle, "Nodes");
        if (nodes != null) {
          for (Element node : elementsOf(nodes, "Object")) {
            nodeItem.put(node.getAttribute("id"), tag);
          }
        }
      }
    }
    System.out.println("item.count=" + items);
    System.out.println("nozzle.count=" + nozzles);

    int segments = 0;
    for (Element object : objects(root)) {
      if (!"Plant/Piping.PipingNetworkSegment".equals(object.getAttribute("type"))) {
        continue;
      }
      segments++;
      String source = nodeItem.getOrDefault(strip(firstReference(object, "SourceNode")), "?");
      String target = nodeItem.getOrDefault(strip(firstReference(object, "TargetNode")), "?");
      System.out.println("segment." + source + "->" + target);
    }
    System.out.println("segment.count=" + segments);

    Dexpi20SemanticValidator.ValidationReport report = Dexpi20SemanticValidator.validate(path);
    System.out.println("semantic.errors=" + report.getErrors().size());
    System.out.println("semantic.warnings=" + report.getWarnings().size());
  }

  private static Iterable<Element> objects(Element root) {
    return elements(root, "Object");
  }

  private static Iterable<Element> elements(Element parent, String tag) {
    NodeList all = parent.getElementsByTagName(tag);
    java.util.List<Element> found = new java.util.ArrayList<>();
    for (int index = 0; index < all.getLength(); index++) {
      found.add((Element) all.item(index));
    }
    return found;
  }

  private static Iterable<Element> elementsOf(Element parent, String tag) {
    NodeList all = parent.getElementsByTagName(tag);
    java.util.List<Element> found = new java.util.ArrayList<>();
    for (int index = 0; index < all.getLength(); index++) {
      org.w3c.dom.Node node = all.item(index);
      if (node.getParentNode() == parent) {
        found.add((Element) node);
      }
    }
    return found;
  }

  private static Element component(Element parent, String property) {
    for (Element element : elementsOf(parent, "Components")) {
      if (property.equals(element.getAttribute("property"))) {
        return element;
      }
    }
    return null;
  }

  private static String data(Element parent, String property) {
    for (Element element : elementsOf(parent, "Data")) {
      if (!property.equals(element.getAttribute("property"))) {
        continue;
      }
      for (Element value : elementsOf(element, "String")) {
        return value.getTextContent().trim();
      }
    }
    return null;
  }

  private static String firstReference(Element parent, String property) {
    for (Element element : elementsOf(parent, "References")) {
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
