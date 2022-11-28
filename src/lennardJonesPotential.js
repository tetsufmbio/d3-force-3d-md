import {binarytree} from "d3-binarytree";
import {quadtree} from "d3-quadtree";
import {octree} from "d3-octree";
import constant from "./constant.js";
import jiggle from "./jiggle.js";
import {x, y, z} from "./simulation.js";
import {atoms} from "./database.js";

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
    treeNode.name = node.name;
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
        const sigma = (atoms[node.name].sigma + atoms[treeNode.name].sigma)/2; // Aqui usei os mesmos dois nodes, mas precisam ser diferentes.
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
    const sigma = (atoms[node.name].sigma + atoms[treeNode.name].sigma)/2; // Aqui usei os mesmos dois nodes, mas precisam ser diferentes.
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