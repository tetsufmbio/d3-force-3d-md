import {binarytree} from "d3-binarytree";
import {quadtree} from "d3-quadtree";
import {octree} from "d3-octree";
import constant from "./constant.js";
import jiggle from "./jiggle.js";
import {x, y, z} from "./simulation.js";

const atoms = {
  ['H']: {mass: 1.00797},
  ['He']: {mass: 4.00260},
  ['Li']: {mass: 6.941},
  ['Be']: {mass: 9.01218},
  ['B']: {mass: 10.81},
  ['C']: {mass: 12.011},
  ['N']: {mass: 14.0067},
  ['O']: {mass: 15.9994},
  ['F']: {mass: 18.998403},
  ['Ne']: {mass: 20.179},
  ['Na']: {mass: 22.98977},
  ['Mg']: {mass: 24.305},
  ['Al']: {mass: 26.98154},
  ['Si']: {mass: 28.0855},
  ['P']: {mass: 30.97376},
  ['S']: {mass: 32.06},
  ['Cl']: {mass: 35.453},
  ['K']: {mass: 39.0983},
  ['Ar']: {mass: 39.948, sigma: 3.405, epsilon: 0.238},
  ['Ca']: {mass: 40.08},
  ['Sc']: {mass: 44.9559},
  ['Ti']: {mass: 47.90},
  ['V']: {mass: 50.9415},
  ['Cr']: {mass: 51.996},
  ['Mn']: {mass: 54.9380},
  ['Fe']: {mass: 55.847},
  ['Ni']: {mass: 58.70},
  ['Co']: {mass: 58.9332},
  ['Cu']: {mass: 63.546},
  ['Zn']: {mass: 65.38},
  ['Ga']: {mass: 69.72},
  ['Ge']: {mass: 72.59},
  ['As']: {mass: 74.9216},
  ['Se']: {mass: 78.96},
  ['Br']: {mass: 79.904},
  ['Kr']: {mass: 83.80},
  ['Rb']: {mass: 85.4678},
  ['Sr']: {mass: 87.62},
  ['Y']: {mass: 88.9059},
  ['Zr']: {mass: 91.22},
  ['Nb']: {mass: 92.9064},
  ['Mo']: {mass: 95.94},
  ['Tc']: {mass: (98)},
  ['Ru']: {mass: 101.07},
  ['Rh']: {mass: 102.9055},
  ['Pd']: {mass: 106.4},
  ['Ag']: {mass: 107.868},
  ['Cd']: {mass: 112.41},
  ['In']: {mass: 114.82},
  ['Sn']: {mass: 118.69},
  ['Sb']: {mass: 121.75},
  ['I']: {mass: 126.9045},
  ['Te']: {mass: 127.60},
  ['Xe']: {mass: 131.30},
  ['Cs']: {mass: 132.9054},
  ['Ba']: {mass: 137.33},
  ['La']: {mass: 138.9055},
  ['Ce']: {mass: 140.12},
  ['Pr']: {mass: 140.9077},
  ['Nd']: {mass: 144.24},
  ['Pm']: {mass: (145)},
  ['Sm']: {mass: 150.4},
  ['Eu']: {mass: 151.96},
  ['Gd']: {mass: 157.25},
  ['Tb']: {mass: 158.9254},
  ['Dy']: {mass: 162.50},
  ['Ho']: {mass: 164.9304},
  ['Er']: {mass: 167.26},
  ['Tm']: {mass: 168.9342},
  ['Yb']: {mass: 173.04},
  ['Lu']: {mass: 174.967},
  ['Hf']: {mass: 178.49},
  ['Ta']: {mass: 180.9479},
  ['W']: {mass: 183.85},
  ['Re']: {mass: 186.207},
  ['Os']: {mass: 190.2},
  ['Ir']: {mass: 192.22},
  ['Pt']: {mass: 195.09},
  ['Au']: {mass: 196.9665},
  ['Hg']: {mass: 200.59},
  ['Tl']: {mass: 204.37},
  ['Pb']: {mass: 207.2},
  ['Bi']: {mass: 208.9804},
  ['Po']: {mass: (209)},
  ['At']: {mass: (210)},
  ['Rn']: {mass: (222)},
  ['Fr']: {mass: (223)},
  ['Ra']: {mass: 226.0254},
  ['Ac']: {mass: 227.0278},
  ['Pa']: {mass: 231.0359},
  ['Th']: {mass: 232.0381},
  ['Np']: {mass: 237.0482},
  ['U']: {mass: 238.029},
  ['Pu']: {mass: (242)},
  ['Am']: {mass: (243)},
  ['Bk']: {mass: (247)},
  ['Cm']: {mass: (247)},
  ['No']: {mass: (250)},
  ['Cf']: {mass: (251)},
  ['Es']: {mass: (252)},
  ['Hs']: {mass: (255)},
  ['Mt']: {mass: (256)},
  ['Fm']: {mass: (257)},
  ['Md']: {mass: (258)},
  ['Lr']: {mass: (260)},
  ['Rf']: {mass: (261)},
  ['Bh']: {mass: (262)},
  ['Db']: {mass: (262)},
  ['Sg']: {mass: (263)},
  ['Uun']: {mass: (269)},
  ['Uuu']: {mass: (272)},
  ['Uub']: {mass: (277)}
}

export default function() {
  var nodes,
      nDim,
      node,
      random,
      alpha,
      strength = constant(1), // strength now serves as a counter
      strengths,
      distanceMin2 = 0.75,
      distanceMax2 = Infinity,
      theta2 = 0.0625, // lennard-jones type interactions mostly short-range
      N=12, M=6; // lennard-jones exponents

  function force(_) {
    var i, 
		n = nodes.length,
		tree = 
			(nDim === 1 ? binarytree(nodes, x)
            :(nDim === 2 ? quadtree(nodes, x, y)
            :(nDim === 3 ? octree(nodes, x, y, z)
            :null
		))).visitAfter(accumulate);
		
    for (alpha = _, i = 0; i < n; ++i) {
      node = nodes[i];
      node.energy = -1; //self-energy
      //node.force_x = 0, node.force_y =0;
      tree.visit(apply);
    }
  }

  function initialize() {
    if (!nodes) return;
    var i, n = nodes.length, node;
    strengths = new Array(n);
    for (i = 0; i < n; ++i) node = nodes[i], strengths[node.index] = +strength(node, i, nodes);
    if (M==N) N += 1; //avoid dividing by zero
  }

  function accumulate(treeNode) {
    var strength = 0, q, c, weight = 0, x, y, z, i;
	var numChildren = treeNode.length;
	
    // For internal nodes, accumulate forces from child quadrants.
    if (numChildren) {
      for (x = y = z = i = 0; i < numChildren; ++i) { // original estah: for (x = y = i = 0; i < 4; ++i)
        if ((q = treeNode[i]) && (c = Math.abs(q.value))) {
          strength += q.value, weight += c, x += c * (q.x || 0), y += c * (q.y || 0), z += c * (q.z || 0);
        }
      }
      strength *= Math.sqrt(4 / numChildren); // scale accumulated strength according to number of dimensions (d3-force-3d)

      treeNode.x = x / weight;
      if (nDim > 1) { treeNode.y = y / weight; }
      if (nDim > 2) { treeNode.z = z / weight; }
    }

    // For leaf nodes, accumulate forces from coincident quadrants.
    else {
      q = treeNode;
      q.x = q.data.x;
      if (nDim > 1) { q.y = q.data.y; }
      if (nDim > 2) { q.z = q.data.z; }
      do strength += strengths[q.data.index];
      while (q = q.next);
    }

    treeNode.value = strength;
  }

  //function apply(treeNode, x1, _, x2) {
  function apply(treeNode, x1, arg1, arg2, arg3) {
    if (!treeNode.value) return true;
    var x2 = [arg1, arg2, arg3][nDim-1];

    var x = treeNode.x - node.x,
        y = (nDim > 1 ? treeNode.y - node.y : 0),
        z = (nDim > 2 ? treeNode.z - node.z : 0),
		w = x2 - x1,
        l = x * x + y * y + z * z,
        force_prefactor;

    // Apply the Barnes-Hut approximation if possible.
    // Limit forces for very close nodes; randomize direction if coincident.
    if (w * w / theta2 < l) {
      if (l < distanceMax2) {
        if (x === 0) x = jiggle(random), l += x * x;
        if (nDim > 1 && y === 0) y = jiggle(random), l += y * y;
        if (nDim > 2 && z === 0) z = jiggle(random), l += z * z;
        if (l < distanceMin2) l = Math.sqrt(distanceMin2 * l);

        const distance = Math.sqrt(l)
        const sigma = (atoms[node.name].sigma + atoms[node.name].sigma)/2; // Aqui usei os mesmos dois nodes, mas precisam ser diferentes.
        const epsilon = 4*atoms[node.name].epsilon;
        force_prefactor = 4*epsilon * (N* Math.pow(sigma/distance, N) - M * Math.pow(sigma/distance, M))/distance;

        node.energy += (N*Math.pow(l,-M/2)-M*Math.pow(l,-N/2))/(M-N)*treeNode.value;
        node.force_x += force_prefactor*x*treeNode.value*alpha;
        if (nDim > 1) { node.force_y += force_prefactor*y*treeNode.value*alpha; }
        if (nDim > 2) { node.force_z += force_prefactor*z*treeNode.value*alpha; }

      }
      return true;
    }

    // Otherwise, process points directly.
    else if (treeNode.length || l >= distanceMax2) return;

    // Limit forces for very close nodes; randomize direction if coincident.
    if (treeNode.data !== node || treeNode.next) {
      if (x === 0) x = jiggle(random), l += x * x;
      if (nDim > 1 && y === 0) y = jiggle(random), l += y * y;
      if (nDim > 2 && z === 0) z = jiggle(random), l += z * z;
      if (l < distanceMin2) l = Math.sqrt(distanceMin2 * l);
    }

    const distance = Math.sqrt(l)
    const sigma = (atoms[node.name].sigma + atoms[node.name].sigma)/2; // Aqui usei os mesmos dois nodes, mas precisam ser diferentes.
    const epsilon = 4*atoms[node.name].epsilon;
    force_prefactor = 4*epsilon * (N* Math.pow(sigma/distance, N) - M * Math.pow(sigma/distance, M))/distance;

    do if (treeNode.data !== node) {
      w = strengths[treeNode.data.index];
      node.energy += (N*Math.pow(l,-M/2)-M*Math.pow(l,-N/2))/(M-N)*w;
      node.force_x += force_prefactor*x*w*alpha;
      if (nDim > 1) { node.force_y += force_prefactor*y*w*alpha; }
      if (nDim > 2) { node.force_z += force_prefactor*z*w*alpha; }

    } while (treeNode = treeNode.next);
  }

  force.initialize = function(_nodes, ...args) {
    nodes = _nodes;
    random = args.find(arg => typeof arg === 'function') || Math.random;
	nDim = args.find(arg => [1, 2, 3].includes(arg)) || 2;
    initialize();
  };

  force.strength = function(_) {
    return arguments.length ? (strength = typeof _ === "function" ? _ : constant(+_), initialize(), force) : strength;
  };

  force.distanceMin = function(_) {
    return arguments.length ? (distanceMin2 = _ * _, force) : Math.sqrt(distanceMin2);
  };

  force.distanceMax = function(_) {
    return arguments.length ? (distanceMax2 = _ * _, force) : Math.sqrt(distanceMax2);
  };

  force.theta = function(_) {
    return arguments.length ? (theta2 = _ * _, force) : Math.sqrt(theta2);
  };
  
  force.repulsivePower = function(_) {
    return arguments.length ? (N = _, force) : N;
  };

  force.attractivePower = function(_) {
    return arguments.length ? (M = _, force) : M;
  };

  return force;
}