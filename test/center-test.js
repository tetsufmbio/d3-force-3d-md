import {forceCenter, forceSimulation} from "../src/index.js";
import {assertNodeEqual} from "./asserts.js";

it("forceCenter repositions nodes", () => {
  const center = forceCenter(0, 0);
  const f = forceSimulation().force("center", center).stop();
  const a = {x: 100, y: 0, z: 0}, b = {x: 200, y: 0, z: 0}, c = {x: 300, y: 0, z: 0};
  f.nodes([a, b, c]);
  f.tick();
  assertNodeEqual(a, {index: 0, x: -100, y: 0, z: 0, vz: 0, vy: 0, vx: 0, force_x: 0, force_y: 0, force_z: 0, mass: 1});
  assertNodeEqual(b, {index: 1, x: 0, y: 0, z: 0, vz: 0, vy: 0, vx: 0, force_x: 0, force_y: 0, force_z: 0, mass: 1});
  assertNodeEqual(c, {index: 2, x: 100, y: 0, z: 0, vz: 0, vy: 0, vx: 0, force_x: 0, force_y: 0, force_z: 0, mass: 1});
});

it("forceCenter respects fixed positions", () => {
  const center = forceCenter();
  const f = forceSimulation().force("center", center).stop();
  const a = {fx: 0, fy: 0, fz: 0}, b = {}, c = {};
  f.nodes([a, b, c]);
  f.tick();
  assertNodeEqual(a, {fx: 0, fy: 0, fz: 0, index: 0, x: 0, y: 0, z: 0, vz: 0, vy: 0, vx: 0, force_x: 0, force_y: 0, force_z: 0, mass: 1});
});
