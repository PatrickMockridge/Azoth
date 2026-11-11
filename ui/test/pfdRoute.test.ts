/**
 * The PFD router's two rules, checked on the path it returns.
 *
 * These are pure geometry tests: the route is a list of orthogonal points, so an obstacle
 * crossing, a diagonal, or a segment retraced in the opposite direction is a fact about the
 * numbers rather than about a screenshot.
 */

import { expect, it } from "vitest";

import { route, type Point, type Rect, type Route } from "../src/pfd/route";

const CLEAR = 12;

function segments(routePath: Route): [Point, Point][] {
  const pairs: [Point, Point][] = [];
  for (let index = 1; index < routePath.points.length; index += 1) {
    const a = routePath.points[index - 1];
    const b = routePath.points[index];
    if (a !== undefined && b !== undefined) {
      pairs.push([a, b]);
    }
  }
  return pairs;
}

function crosses(a: Point, b: Point, rect: Rect): boolean {
  const left = rect.x;
  const right = rect.x + rect.width;
  const top = rect.y;
  const bottom = rect.y + rect.height;
  if (a.y === b.y && a.y > top && a.y < bottom) {
    return Math.max(a.x, b.x) > left && Math.min(a.x, b.x) < right;
  }
  if (a.x === b.x && a.x > left && a.x < right) {
    return Math.max(a.y, b.y) > top && Math.min(a.y, b.y) < bottom;
  }
  return false;
}

function overlaps(a: [Point, Point], b: [Point, Point]): boolean {
  const [a1, a2] = a;
  const [b1, b2] = b;
  if (a1.y === a2.y && b1.y === b2.y && a1.y === b1.y) {
    return Math.max(Math.min(a1.x, a2.x), Math.min(b1.x, b2.x)) <
      Math.min(Math.max(a1.x, a2.x), Math.max(b1.x, b2.x));
  }
  if (a1.x === a2.x && b1.x === b2.x && a1.x === b1.x) {
    return Math.max(Math.min(a1.y, a2.y), Math.min(b1.y, b2.y)) <
      Math.min(Math.max(a1.y, a2.y), Math.max(b1.y, b2.y));
  }
  return false;
}

it("routes a forward edge around a plate without touching it", () => {
  const obstacle = { x: 150, y: 60, width: 100, height: 80 };
  const path = route({
    sx: 0,
    sy: 100,
    tx: 400,
    ty: 100,
    sourceSide: "right",
    targetSide: "left",
    obstacles: [obstacle],
    recycle: false,
    lane: 0,
  });
  expect(path).not.toBeNull();
  if (path === null) {
    return;
  }
  for (const [a, b] of segments(path)) {
    expect(a.x === b.x || a.y === b.y, "a segment is diagonal").toBe(true);
    expect(crosses(a, b, obstacle), "a segment crosses the plate").toBe(false);
  }
});

it("returns a recycle below the plates, not through them", () => {
  const plates: Rect[] = [
    { x: 0, y: 100, width: 100, height: 100 },
    { x: 150, y: 100, width: 100, height: 100 },
    { x: 300, y: 100, width: 100, height: 100 },
  ];
  const path = route({
    sx: 400,
    sy: 200,
    tx: 0,
    ty: 150,
    sourceSide: "bottom",
    targetSide: "left",
    obstacles: plates,
    recycle: true,
    lane: 0,
  });
  expect(path).not.toBeNull();
  if (path === null) {
    return;
  }
  for (const [a, b] of segments(path)) {
    expect(a.x === b.x || a.y === b.y, "a segment is diagonal").toBe(true);
    for (const plate of plates) {
      expect(crosses(a, b, plate), "a segment crosses a plate").toBe(false);
    }
  }
  const bottom = Math.max(...plates.map((plate) => plate.y + plate.height));
  const lane = segments(path).some(
    ([a, b]) => a.y === b.y && a.y > bottom + CLEAR && Math.abs(a.x - b.x) > 0,
  );
  expect(lane, "the return does not use a lane below the plates").toBe(true);
});

it("keeps two recycles on different return lanes", () => {
  const plates: Rect[] = [
    { x: 0, y: 100, width: 100, height: 100 },
    { x: 250, y: 100, width: 100, height: 100 },
  ];
  const first = route({
    sx: 350,
    sy: 200,
    tx: 0,
    ty: 150,
    sourceSide: "bottom",
    targetSide: "left",
    obstacles: plates,
    recycle: true,
    lane: 0,
  });
  const second = route({
    sx: 350,
    sy: 200,
    tx: 0,
    ty: 150,
    sourceSide: "bottom",
    targetSide: "left",
    obstacles: plates,
    recycle: true,
    lane: 1,
  });
  expect(first).not.toBeNull();
  expect(second).not.toBeNull();
  if (first === null || second === null) {
    return;
  }
  const laneOf = (routePath: Route): number | null => {
    const bottom = Math.max(...plates.map((plate) => plate.y + plate.height));
    const horizontal = segments(routePath).find(
      ([a, b]) => a.y === b.y && a.y > bottom + CLEAR && Math.abs(a.x - b.x) > 0,
    );
    return horizontal?.[0].y ?? null;
  };
  expect(laneOf(first)).not.toBe(laneOf(second));
});

it("does not retrace a line in the opposite direction", () => {
  const plates: Rect[] = [
    { x: 0, y: 100, width: 100, height: 100 },
    { x: 250, y: 100, width: 100, height: 100 },
  ];
  const path = route({
    sx: 350,
    sy: 200,
    tx: 0,
    ty: 150,
    sourceSide: "bottom",
    targetSide: "left",
    obstacles: plates,
    recycle: true,
    lane: 0,
  });
  expect(path).not.toBeNull();
  if (path === null) {
    return;
  }
  const all = segments(path);
  for (let index = 0; index < all.length; index += 1) {
    for (let other = index + 2; other < all.length; other += 1) {
      expect(overlaps(all[index] as [Point, Point], all[other] as [Point, Point])).toBe(
        false,
      );
    }
  }
});
