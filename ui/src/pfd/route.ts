/**
 * Orthogonal edge routing around the node plates.
 *
 * **Two rules, checked rather than hoped for:**
 *
 * 1. A connector never crosses a node plate. Every plate is an obstacle, expanded by a clearance;
 *    the router only accepts segments that stay outside all of them.
 * 2. A connector never runs along itself in two directions. The route is a shortest path on a
 *    visibility grid with a turn cost, and a recycle is pushed to a dedicated lane below the
 *    drawing - so the forward and return halves occupy different coordinates.
 *
 * **A small grid, not a full geometry library.** The candidates are the plate edges, the midpoints
 * between them, the two stubs and the lane; Dijkstra over their product is a few hundred states.
 * The path is then simplified to its bends, which is the SVG `path` and the same list the tests
 * check for crossings and reversals.
 */

export type Side = "left" | "right" | "top" | "bottom";

export interface Point {
  x: number;
  y: number;
}

export interface Rect {
  x: number;
  y: number;
  width: number;
  height: number;
}

export interface RouteInput {
  sx: number;
  sy: number;
  tx: number;
  ty: number;
  sourceSide: Side;
  targetSide: Side;
  obstacles: readonly Rect[];
  recycle: boolean;
  /** The index of the return lane; only meaningful for a recycle. */
  lane: number;
}

export interface Route {
  points: Point[];
  path: string;
  /** The point halfway along the route, for the label. */
  mid: Point;
}

/** Clearance between a plate and a line that passes it. */
const CLEAR = 12;
/** How far a line leaves a handle before it is allowed to turn. */
const STUB = 18;
/** The first return lane's distance below the lowest plate, and the spacing between lanes. */
const LANE_GAP = 28;
const LANE_STEP = 18;
/** A bend costs this much in the search, so the result is a few long runs rather than a staircase. */
const TURN = 40;

function step(side: Side): Point {
  switch (side) {
    case "left":
      return { x: -1, y: 0 };
    case "right":
      return { x: 1, y: 0 };
    case "top":
      return { x: 0, y: -1 };
    case "bottom":
      return { x: 0, y: 1 };
  }
}

function blocked(point: Point, obstacles: readonly Rect[]): boolean {
  return obstacles.some(
    (rect) =>
      point.x > rect.x - CLEAR &&
      point.x < rect.x + rect.width + CLEAR &&
      point.y > rect.y - CLEAR &&
      point.y < rect.y + rect.height + CLEAR,
  );
}

/** Whether an axis-aligned segment passes through a plate's clearance box. */
function crosses(a: Point, b: Point, obstacles: readonly Rect[]): boolean {
  for (const rect of obstacles) {
    const left = rect.x - CLEAR;
    const right = rect.x + rect.width + CLEAR;
    const top = rect.y - CLEAR;
    const bottom = rect.y + rect.height + CLEAR;
    if (a.y === b.y && a.y > top && a.y < bottom) {
      if (Math.max(a.x, b.x) > left && Math.min(a.x, b.x) < right) {
        return true;
      }
    }
    if (a.x === b.x && a.x > left && a.x < right) {
      if (Math.max(a.y, b.y) > top && Math.min(a.y, b.y) < bottom) {
        return true;
      }
    }
  }
  return false;
}

function sorted(values: Set<number>): number[] {
  return [...values].sort((a, b) => a - b);
}

/** Add the midpoint of every adjacent pair, which is where a channel between two plates lives. */
function withMidpoints(values: number[]): number[] {
  const result = [...values];
  for (let index = 1; index < values.length; index += 1) {
    const a = values[index - 1] ?? 0;
    const b = values[index] ?? 0;
    result.push((a + b) / 2);
  }
  return sorted(new Set(result));
}

function simplify(points: Point[]): Point[] {
  const unique: Point[] = [];
  for (const point of points) {
    const last = unique.at(-1);
    if (last === undefined || last.x !== point.x || last.y !== point.y) {
      unique.push(point);
    }
  }
  const result: Point[] = [];
  for (const point of unique) {
    while (result.length >= 2) {
      const a = result[result.length - 2];
      const b = result[result.length - 1];
      if (a === undefined || b === undefined) {
        break;
      }
      const collinear = (a.x === b.x && b.x === point.x) || (a.y === b.y && b.y === point.y);
      if (!collinear) {
        break;
      }
      result.pop();
    }
    result.push(point);
  }
  return result;
}

function asPath(points: Point[]): string {
  return points
    .map((point, index) => `${index === 0 ? "M" : "L"} ${point.x} ${point.y}`)
    .join(" ");
}

function midpoint(points: Point[]): Point {
  let total = 0;
  for (let index = 1; index < points.length; index += 1) {
    const a = points[index - 1];
    const b = points[index];
    if (a !== undefined && b !== undefined) {
      total += Math.abs(a.x - b.x) + Math.abs(a.y - b.y);
    }
  }
  let walked = 0;
  for (let index = 1; index < points.length; index += 1) {
    const a = points[index - 1];
    const b = points[index];
    if (a === undefined || b === undefined) {
      continue;
    }
    const length = Math.abs(a.x - b.x) + Math.abs(a.y - b.y);
    if (walked + length >= total / 2) {
      const along = total === 0 ? 0 : (total / 2 - walked) / length;
      return { x: a.x + (b.x - a.x) * along, y: a.y + (b.y - a.y) * along };
    }
    walked += length;
  }
  return points[0] ?? { x: 0, y: 0 };
}

/**
 * Route one edge.
 *
 * Returns `null` when the grid offers no path - a stack of plates with no gap is a drawing no
 * orthogonal line can honour - and the caller falls back to xyflow's smooth step for that edge.
 */
export function route(input: RouteInput): Route | null {
  const out = step(input.sourceSide);
  const into = step(input.targetSide);
  const start = { x: input.sx + out.x * STUB, y: input.sy + out.y * STUB };
  const end = { x: input.tx + into.x * STUB, y: input.ty + into.y * STUB };

  if (input.obstacles.length === 0) {
    const points = simplify([
      { x: input.sx, y: input.sy },
      start,
      end,
      { x: input.tx, y: input.ty },
    ]);
    return { points, path: asPath(points), mid: midpoint(points) };
  }

  const xs = new Set<number>([input.sx, input.tx, start.x, end.x]);
  const ys = new Set<number>([input.sy, input.ty, start.y, end.y]);
  for (const rect of input.obstacles) {
    xs.add(rect.x - CLEAR);
    xs.add(rect.x + rect.width + CLEAR);
    ys.add(rect.y - CLEAR);
    ys.add(rect.y + rect.height + CLEAR);
  }
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);
  xs.add(minX - STUB);
  xs.add(maxX + STUB);
  ys.add(minY - STUB);
  ys.add(maxY + STUB);

  // The return lanes are below every plate, one per recycle, so they cannot share a line.
  const laneY = maxY + LANE_GAP + Math.max(0, input.lane) * LANE_STEP;
  if (input.recycle) {
    ys.add(laneY);
  }

  const xList = withMidpoints(sorted(xs));
  const yList = withMidpoints(sorted(ys));
  const xi = new Map(xList.map((value, index) => [value, index]));
  const yi = new Map(yList.map((value, index) => [value, index]));

  const startXi = xi.get(start.x);
  const startYi = yi.get(start.y);
  const endXi = xi.get(end.x);
  const endYi = yi.get(end.y);
  if (
    startXi === undefined ||
    startYi === undefined ||
    endXi === undefined ||
    endYi === undefined ||
    blocked(start, input.obstacles) ||
    blocked(end, input.obstacles)
  ) {
    return null;
  }

  const width = xList.length;
  const height = yList.length;
  const at = (x: number, y: number): Point => ({
    x: xList[x] ?? 0,
    y: yList[y] ?? 0,
  });

  type State = { x: number; y: number; dir: number };
  const key = (state: State): string => `${state.x},${state.y},${state.dir}`;
  const distance = new Map<string, number>();
  const previous = new Map<string, State>();
  type Entry = { state: State; cost: number; priority: number };
  const heap: Entry[] = [];
  const earlier = (a: Entry, b: Entry): boolean => a.priority < b.priority;
  const push = (entry: Entry): void => {
    heap.push(entry);
    let index = heap.length - 1;
    while (index > 0) {
      const parent = Math.floor((index - 1) / 2);
      const parentEntry = heap[parent];
      const entryHere = heap[index];
      if (parentEntry === undefined || entryHere === undefined || !earlier(entryHere, parentEntry)) {
        break;
      }
      heap[parent] = entryHere;
      heap[index] = parentEntry;
      index = parent;
    }
  };
  const pop = (): Entry | undefined => {
    const top = heap[0];
    const last = heap.pop();
    if (top === undefined || last === undefined) {
      return top;
    }
    if (heap.length > 0) {
      heap[0] = last;
      let index = 0;
      while (true) {
        const left = index * 2 + 1;
        const right = left + 1;
        let smallest = index;
        for (const candidate of [left, right]) {
          const here = heap[candidate];
          const best = heap[smallest];
          if (here !== undefined && best !== undefined && earlier(here, best)) {
            smallest = candidate;
          }
        }
        if (smallest === index) {
          break;
        }
        const a = heap[index];
        const b = heap[smallest];
        if (a === undefined || b === undefined) {
          break;
        }
        heap[index] = b;
        heap[smallest] = a;
        index = smallest;
      }
    }
    return top;
  };
  const heuristic = (point: Point): number =>
    Math.abs(point.x - (xList[endXi] ?? 0)) + Math.abs(point.y - (yList[endYi] ?? 0));
  const startState: State = { x: startXi, y: startYi, dir: -1 };
  distance.set(key(startState), 0);
  push({ state: startState, cost: 0, priority: heuristic(start) });

  const moves = [
    { dx: 1, dy: 0, dir: 0 },
    { dx: -1, dy: 0, dir: 1 },
    { dx: 0, dy: 1, dir: 2 },
    { dx: 0, dy: -1, dir: 3 },
  ];
  let goal: State | null = null;

  while (heap.length > 0) {
    const next = pop();
    if (next === undefined) {
      break;
    }
    const current = next.state;
    if (next.cost > (distance.get(key(current)) ?? Number.POSITIVE_INFINITY)) {
      continue;
    }
    if (current.x === endXi && current.y === endYi) {
      goal = current;
      break;
    }
    const a = at(current.x, current.y);
    for (const move of moves) {
      const x = current.x + move.dx;
      const y = current.y + move.dy;
      if (x < 0 || y < 0 || x >= width || y >= height) {
        continue;
      }
      const b = at(x, y);
      if (crosses(a, b, input.obstacles)) {
        continue;
      }
      // A recycle's horizontal runs belong in its return lane; moving sideways in the node band
      // is expensive rather than forbidden, so a pathological layout still gets a path.
      const sideways = move.dy === 0;
      const aboveLane = b.y < laneY - 0.01;
      const penalty = input.recycle && sideways && aboveLane ? 400 : 0;
      const turn = current.dir !== -1 && current.dir !== move.dir ? TURN : 0;
      const cost =
        next.cost + Math.abs(a.x - b.x) + Math.abs(a.y - b.y) + turn + penalty;
      const state: State = { x, y, dir: move.dir };
      if (cost < (distance.get(key(state)) ?? Number.POSITIVE_INFINITY)) {
        distance.set(key(state), cost);
        previous.set(key(state), current);
        push({ state, cost, priority: cost + heuristic(b) });
      }
    }
  }

  if (goal === null) {
    return null;
  }

  const grid: Point[] = [];
  let cursor: State | undefined = goal;
  while (cursor !== undefined) {
    grid.push(at(cursor.x, cursor.y));
    cursor = previous.get(key(cursor));
  }
  grid.reverse();

  const points = simplify([
    { x: input.sx, y: input.sy },
    ...grid,
    { x: input.tx, y: input.ty },
  ]);
  return { points, path: asPath(points), mid: midpoint(points) };
}
