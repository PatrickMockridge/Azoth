/**
 * The worked cases the editor offers.
 *
 * **The documents are imported, not copied.** A copy is a second flowsheet that can go stale
 * against the one the Rust tests run; the glob makes every new file in the directory visible to
 * the menu in the same build that ships it. The descriptions are the one hand-written part, and
 * the test beside this module holds them to the directory one for one.
 */

const documents = import.meta.glob<string>("../../specs/flowsheets/examples/*.toml", {
  eager: true,
  query: "?raw",
  import: "default",
});

/** A case the welcome card and the Case menu can open. */
export interface Example {
  /** The file name, which is stable and is the menu's value. */
  file: string;
  /** Its title, as the menu shows it. */
  title: string;
  /** One sentence on what it demonstrates. */
  description: string;
  /** The document itself. */
  document: string;
}

/** What each shipped example is for. Every file needs an entry; the test says so. */
const META: Record<string, { title: string; description: string }> = {
  "heat_a_gas.toml": {
    title: "Heat a gas",
    description: "One feed, one heater and one product: the smallest complete case.",
  },
  "two_feed_mixer.toml": {
    title: "Two feeds, mix and recycle",
    description:
      "Two distinct feeds meet at a mixer with a separate inlet for each, then a pump, heater, separator and live liquid recycle.",
  },
  "compress_and_cool.toml": {
    title: "Compress, cool and separate",
    description: "A rich gas is compressed, cooled and split into sales gas and condensate.",
  },
  "purge_and_recycle.toml": {
    title: "Split, purge and recycle",
    description: "A splitter closes a recycle: half the flow is purged, half returns.",
  },
  "pipe_and_valve.toml": {
    title: "Pipe, valve and cold separation",
    description: "A pipeline drops pressure through a valve and condenses in a cold separator.",
  },
};

/** The examples, in the order the menu lists them. */
export const EXAMPLES: Example[] = Object.entries(documents)
  .map(([path, document]) => {
    const file = path.split("/").pop() ?? path;
    const meta = META[file];
    return {
      file,
      // A file without a description is a menu entry nobody thought about, which the test
      // refuses; this fallback keeps the page drawable while that test says so.
      title: meta?.title ?? file.replace(/\.toml$/, "").replace(/_/g, " "),
      description: meta?.description ?? "",
      document,
    };
  })
  .sort((left, right) => left.title.localeCompare(right.title));
