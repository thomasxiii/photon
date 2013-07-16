/*
 * Photon
 * http://photon.attasi.com
 *
 * Licensed under the MIT license.
 * Copyright 2012 Tom Giannattasio
 */



var Photon = {
	version: '0.0.3',

	degToRad: function(deg) {
		return deg * Math.PI / 180;
	},

	radToDeg: function(rad) {
		return rad * 180 / Math.PI;
	},

	getRotationVector: function(originVector, rotations) {
		var xVector = originVector.rotate(rotations.x, Line.create([0, 0, 0], [1, 0, 0]));
		var yVector = xVector.rotate(rotations.y, Line.create([0, 0, 0], [0, 1, 0]));
		var zVector = yVector.rotate(rotations.z, Line.create([0, 0, 0], [0, 0, 1]));
		return zVector;
	},

  getTransformString: function() {
    if(Photon.transformString) {
      return Photon.transformString;
    }

    var transformString;
    var tests = ['transform', 'webkitTransform', 'MozTransform', 'msTransform', 'OTransform'];
    var element = document.createElement('div');

    for(var i = 0; i < tests.length; i++) {
      if(element.style[tests[i]] == '') {
        transformString = tests[i];
      }
    }

    Photon.transformString = transformString;
    return transformString;
  },

	// converts transform matrix into a WebKitCSSMatrix object.
	// multiplies values to avoid whackification
	buildMatrix: function(faceTransform) {
		var matrix = new FirminCSSMatrix(faceTransform);
		
		matrix.m11 = matrix.m11 * 1e16;
		matrix.m12 = matrix.m12 * 1e16;
		matrix.m13 = matrix.m13 * 1e16;
		matrix.m14 = matrix.m14 * 1e16;
		
		matrix.m21 = matrix.m21 * 1e16;
		matrix.m22 = matrix.m22 * 1e16;
		matrix.m23 = matrix.m23 * 1e16;
		matrix.m24 = matrix.m24 * 1e16;
		
		matrix.m31 = matrix.m31 * 1e16;
		matrix.m32 = matrix.m32 * 1e16;
		matrix.m33 = matrix.m33 * 1e16;
		matrix.m34 = matrix.m34 * 1e16;
		
		matrix.m41 = matrix.m41 * 1e16;
		matrix.m42 = matrix.m42 * 1e16;
		matrix.m43 = matrix.m43 * 1e16;
		matrix.m44 = matrix.m44 * 1e16;
		
		return matrix;
	}
};



Photon.Light = function(xVal, yVal, zVal) {
	this.moveTo(xVal || 0, yVal || 0, zVal || 100);
	this.calculateVector();
}

Photon.Light.prototype = {
	moveTo: function(x, y, z) {
		this.x = x;
		this.y = y;
		this.z = z;
		this.calculateVector();
	},

	// covert the light coordinates into a vector
	calculateVector: function() {
		this.magnitude = Math.sqrt((this.x * this.x) + (this.y * this.y) + (this.z * this.z));
		this.vector = $V([this.x / this.magnitude, this.y / this.magnitude, this.z / this.magnitude]);
	}
}



Photon.Face = function(element, maxShade, maxTint, isBackfaced) {
	// set properties
	this.element = element;
	this.maxShade = maxShade || .5;
	this.maxTint = maxTint || 0;
	this.isBackfaced = isBackfaced || false;

	// create shader element
	this.shaderElement = new Photon.ShaderElement(this.element);
	this.element.insertBefore(this.shaderElement, this.element.firstChild);
	
  this.transformString = Photon.getTransformString();

	// calculate absolute rotations
	this.getRotations();
}

Photon.Face.prototype = {
	getRotations: function() {
		// pull the transform property
    var faceTransform = window.getComputedStyle(this.element)[this.transformString] || 'matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)';

		// convert the transform data into a matrix
		this.matrix = Photon.buildMatrix(faceTransform);

		// extract the individual transform values
		var faceDecomp = this.matrix.decompose();
		this.rotations = {
			x: faceDecomp.rotate.x,
			y: faceDecomp.rotate.y,
			z: faceDecomp.rotate.z
		};
		
		// set the face's vector
		this.vector = Photon.getRotationVector($V([0, 0, 1]), this.rotations);
	},
	
	render: function(light, getNewRotations, parentRotations) {
		if(getNewRotations) {
			this.getRotations();
		}
		
		// calculate the absolute vector
		var fullVector;
		if(parentRotations) {
			fullVector = Photon.getRotationVector(this.vector, parentRotations);
		} else {
			fullVector = this.vector;
		}

		// calculate the anglar distance from the light
		this.angleFrom = Photon.radToDeg(light.vector.angleFrom(fullVector));

		// determine the background color of the shader element
		var background;
		
		// var anglePercentage = this.angleFrom / 180;
		var anglePercentage = this.isBackfaced ? this.angleFrom / 180 : this.angleFrom / 90;
		if(this.isBackfaced && anglePercentage > .5) {
			anglePercentage = 1 - anglePercentage;
		}
		var range = Math.abs(this.maxShade + this.maxTint);
		var rangedPercentage = range * anglePercentage;
		this.rangedPercentage = rangedPercentage;

		// determine whether to shade or tint
		if(rangedPercentage <= this.maxTint) {
			background = 'rgba(255, 255, 255, ' + Math.abs(this.maxTint - rangedPercentage) + ')';
		} else {
			background = 'rgba(0, 0, 0, ' + Math.abs(rangedPercentage - this.maxTint) + ')';
		}

		// apply the shading
		this.shaderElement.style.background = background;
	},

	setMaxShade: function(value) {
		this.maxShade = value;
	},

	setMaxTint: function(value) {
		this.maxTint = value;
	}
};


// create the element to used for shading and tinting
Photon.ShaderElement = function(parent) {
	var shaderElement = document.createElement('div');
	shaderElement.className = 'photon-shader';
	shaderElement.style.position = 'absolute';
	shaderElement.style.top = '0';
	shaderElement.style.left = '0';
	shaderElement.style.width = window.getComputedStyle(parent).width;
	shaderElement.style.height = window.getComputedStyle(parent).height;

	return shaderElement;
}


// a group of faces within a single parent object
Photon.FaceGroup = function(parent, faces, maxShade, maxTint, isBackfaced) {
	this.element = parent;
	this.faces = [];
  this.transformString = Photon.getTransformString();

	var childFaces = faces;
	for(var i = 0; i < childFaces.length; i++) {
		this.faces[i] = new Photon.Face(childFaces[i], maxShade, maxTint, isBackfaced);
	}
}

Photon.FaceGroup.prototype = {
	getRotations: function() {
		var faceTransform = window.getComputedStyle(this.element)[this.transformString] || 'matrix3d(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)';

		this.matrix = Photon.buildMatrix(faceTransform);
		var faceDecomp = this.matrix.decompose();
		
		this.rotations = {
			x: faceDecomp.rotate.x,
			y: faceDecomp.rotate.y,
			z: faceDecomp.rotate.z
		};
		
		this.vector = Photon.getRotationVector($V([0, 0, 1]), this.rotations);
	},

	render: function(light, getNewGroupRotations, getNewFaceRotations) {
		if(getNewGroupRotations) {
			this.getRotations();
		}
		
		this.angleFrom = Photon.radToDeg(light.vector.angleFrom(this.vector));

		for(var i = 0, length = this.faces.length; i < length; i++) {
			this.faces[i].render(light, getNewFaceRotations, this.rotations);
		}
	},

	setMaxShade: function(value) {
		for(var i = 0; i < this.faces.length; i++) {
			this.faces[i].setMaxShade(value);
		}
	},

	setMaxTint: function(value) {
		for(var i = 0; i < this.faces.length; i++) {
			this.faces[i].setMaxTint(value);
		}
	}
};












// === Sylvester ===
// Vector and Matrix mathematics modules for JavaScript
// Copyright (c) 2007 James Coglan
// 
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.


var Sylvester = {
  version: '0.1.3',
  precision: 1e-6
};



function Vector() {}
Vector.prototype = {

  // Returns the modulus ('length') of the vector
  modulus: function() {
    return Math.sqrt(this.dot(this));
  },

  // Returns a copy of the vector
  dup: function() {
    return Vector.create(this.elements);
  },
  
  // Calls the iterator for each element of the vector in turn
  each: function(fn) {
    var n = this.elements.length, k = n, i;
    do { i = k - n;
      fn(this.elements[i], i+1);
    } while (--n);
  },

  // Returns the angle between the vector and the argument (also a vector)
  angleFrom: function(vector) {
    var V = vector.elements || vector;
    var n = this.elements.length, k = n, i;
    if (n != V.length) { return null; }
    var dot = 0, mod1 = 0, mod2 = 0;
    // Work things out in parallel to save time
    this.each(function(x, i) {
      dot += x * V[i-1];
      mod1 += x * x;
      mod2 += V[i-1] * V[i-1];
    });
    mod1 = Math.sqrt(mod1); mod2 = Math.sqrt(mod2);
    if (mod1*mod2 === 0) { return null; }
    var theta = dot / (mod1*mod2);
    if (theta < -1) { theta = -1; }
    if (theta > 1) { theta = 1; }
    return Math.acos(theta);
  },

  // Returns the scalar product of the vector with the argument
  // Both vectors must have equal dimensionality
  dot: function(vector) {
    var V = vector.elements || vector;
    var i, product = 0, n = this.elements.length;
    if (n != V.length) { return null; }
    do { product += this.elements[n-1] * V[n-1]; } while (--n);
    return product;
  },

  // Rotates the vector about the given object. The object should be a 
  // point if the vector is 2D, and a line if it is 3D. Be careful with line directions!
  rotate: function(t, obj) {
    var V, R, x, y, z;
    switch (this.elements.length) {
      case 2:
        V = obj.elements || obj;
        if (V.length != 2) { return null; }
        R = Matrix.Rotation(t).elements;
        x = this.elements[0] - V[0];
        y = this.elements[1] - V[1];
        return Vector.create([
          V[0] + R[0][0] * x + R[0][1] * y,
          V[1] + R[1][0] * x + R[1][1] * y
        ]);
        break;
      case 3:
        if (!obj.direction) { return null; }
        var C = obj.pointClosestTo(this).elements;
        R = Matrix.Rotation(t, obj.direction).elements;
        x = this.elements[0] - C[0];
        y = this.elements[1] - C[1];
        z = this.elements[2] - C[2];
        return Vector.create([
          C[0] + R[0][0] * x + R[0][1] * y + R[0][2] * z,
          C[1] + R[1][0] * x + R[1][1] * y + R[1][2] * z,
          C[2] + R[2][0] * x + R[2][1] * y + R[2][2] * z
        ]);
        break;
      default:
        return null;
    }
  },

  // Set vector's elements from an array
  setElements: function(els) {
    this.elements = (els.elements || els).slice();
    return this;
  }
};

// Constructor function
Vector.create = function(elements) {
  var V = new Vector();
  return V.setElements(elements);
};

var $V = Vector.create;



function Line() {}
Line.prototype = {

  // Returns the line's perpendicular distance from the argument,
  // which can be a point, a line or a plane
  distanceFrom: function(obj) {
    if (obj.normal) { return obj.distanceFrom(this); }
    if (obj.direction) {
      // obj is a line
      if (this.isParallelTo(obj)) { return this.distanceFrom(obj.anchor); }
      var N = this.direction.cross(obj.direction).toUnitVector().elements;
      var A = this.anchor.elements, B = obj.anchor.elements;
      return Math.abs((A[0] - B[0]) * N[0] + (A[1] - B[1]) * N[1] + (A[2] - B[2]) * N[2]);
    } else {
      // obj is a point
      var P = obj.elements || obj;
      var A = this.anchor.elements, D = this.direction.elements;
      var PA1 = P[0] - A[0], PA2 = P[1] - A[1], PA3 = (P[2] || 0) - A[2];
      var modPA = Math.sqrt(PA1*PA1 + PA2*PA2 + PA3*PA3);
      if (modPA === 0) return 0;
      // Assumes direction vector is normalized
      var cosTheta = (PA1 * D[0] + PA2 * D[1] + PA3 * D[2]) / modPA;
      var sin2 = 1 - cosTheta*cosTheta;
      return Math.abs(modPA * Math.sqrt(sin2 < 0 ? 0 : sin2));
    }
  },

  // Returns true iff the argument is a point on the line
  contains: function(point) {
    var dist = this.distanceFrom(point);
    return (dist !== null && dist <= Sylvester.precision);
  },

  // Returns the point on the line that is closest to the given point or line
  pointClosestTo: function(obj) {
    if (obj.direction) {
      // obj is a line
      if (this.intersects(obj)) { return this.intersectionWith(obj); }
      if (this.isParallelTo(obj)) { return null; }
      var D = this.direction.elements, E = obj.direction.elements;
      var D1 = D[0], D2 = D[1], D3 = D[2], E1 = E[0], E2 = E[1], E3 = E[2];
      // Create plane containing obj and the shared normal and intersect this with it
      // Thank you: http://www.cgafaq.info/wiki/Line-line_distance
      var x = (D3 * E1 - D1 * E3), y = (D1 * E2 - D2 * E1), z = (D2 * E3 - D3 * E2);
      var N = Vector.create([x * E3 - y * E2, y * E1 - z * E3, z * E2 - x * E1]);
      var P = Plane.create(obj.anchor, N);
      return P.intersectionWith(this);
    } else {
      // obj is a point
      var P = obj.elements || obj;
      if (this.contains(P)) { return Vector.create(P); }
      var A = this.anchor.elements, D = this.direction.elements;
      var D1 = D[0], D2 = D[1], D3 = D[2], A1 = A[0], A2 = A[1], A3 = A[2];
      var x = D1 * (P[1]-A2) - D2 * (P[0]-A1), y = D2 * ((P[2] || 0) - A3) - D3 * (P[1]-A2),
          z = D3 * (P[0]-A1) - D1 * ((P[2] || 0) - A3);
      var V = Vector.create([D2 * x - D3 * z, D3 * y - D1 * x, D1 * z - D2 * y]);
      var k = this.distanceFrom(P) / V.modulus();
      return Vector.create([
        P[0] + V.elements[0] * k,
        P[1] + V.elements[1] * k,
        (P[2] || 0) + V.elements[2] * k
      ]);
    }
  },

  // Returns a copy of the line rotated by t radians about the given line. Works by
  // finding the argument's closest point to this line's anchor point (call this C) and
  // rotating the anchor about C. Also rotates the line's direction about the argument's.
  // Be careful with this - the rotation axis' direction affects the outcome!
  rotate: function(t, line) {
    // If we're working in 2D
    if (typeof(line.direction) == 'undefined') { line = Line.create(line.to3D(), Vector.k); }
    var R = Matrix.Rotation(t, line.direction).elements;
    var C = line.pointClosestTo(this.anchor).elements;
    var A = this.anchor.elements, D = this.direction.elements;
    var C1 = C[0], C2 = C[1], C3 = C[2], A1 = A[0], A2 = A[1], A3 = A[2];
    var x = A1 - C1, y = A2 - C2, z = A3 - C3;
    return Line.create([
      C1 + R[0][0] * x + R[0][1] * y + R[0][2] * z,
      C2 + R[1][0] * x + R[1][1] * y + R[1][2] * z,
      C3 + R[2][0] * x + R[2][1] * y + R[2][2] * z
    ], [
      R[0][0] * D[0] + R[0][1] * D[1] + R[0][2] * D[2],
      R[1][0] * D[0] + R[1][1] * D[1] + R[1][2] * D[2],
      R[2][0] * D[0] + R[2][1] * D[1] + R[2][2] * D[2]
    ]);
  },

  // Set the line's anchor point and direction.
  setVectors: function(anchor, direction) {
    // Need to do this so that line's properties are not
    // references to the arguments passed in
    anchor = Vector.create(anchor);
    direction = Vector.create(direction);
    if (anchor.elements.length == 2) {anchor.elements.push(0); }
    if (direction.elements.length == 2) { direction.elements.push(0); }
    if (anchor.elements.length > 3 || direction.elements.length > 3) { return null; }
    var mod = direction.modulus();
    if (mod === 0) { return null; }
    this.anchor = anchor;
    this.direction = Vector.create([
      direction.elements[0] / mod,
      direction.elements[1] / mod,
      direction.elements[2] / mod
    ]);
    return this;
  }
};

// Constructor function
Line.create = function(anchor, direction) {
  var L = new Line();
  return L.setVectors(anchor, direction);
};



function Matrix() {}
Matrix.prototype = {
  // Set the matrix's elements from an array. If the argument passed
  // is a vector, the resulting matrix will be a single column.
  setElements: function(els) {
    var i, elements = els.elements || els;
    if (typeof(elements[0][0]) != 'undefined') {
      var ni = elements.length, ki = ni, nj, kj, j;
      this.elements = [];
      do { i = ki - ni;
        nj = elements[i].length; kj = nj;
        this.elements[i] = [];
        do { j = kj - nj;
          this.elements[i][j] = elements[i][j];
        } while (--nj);
      } while(--ni);
      return this;
    }
    var n = elements.length, k = n;
    this.elements = [];
    do { i = k - n;
      this.elements.push([elements[i]]);
    } while (--n);
    return this;
  }
};

// Constructor function
Matrix.create = function(elements) {
  var M = new Matrix();
  return M.setElements(elements);
};

Matrix.Rotation = function(theta, a) {
  if (!a) {
    return Matrix.create([
      [Math.cos(theta),  -Math.sin(theta)],
      [Math.sin(theta),   Math.cos(theta)]
    ]);
  }
  var axis = a.dup();
  if (axis.elements.length != 3) { return null; }
  var mod = axis.modulus();
  var x = axis.elements[0]/mod, y = axis.elements[1]/mod, z = axis.elements[2]/mod;
  var s = Math.sin(theta), c = Math.cos(theta), t = 1 - c;
  // Formula derived here: http://www.gamedev.net/reference/articles/article1199.asp
  // That proof rotates the co-ordinate system so theta
  // becomes -theta and sin becomes -sin here.
  return Matrix.create([
    [ t*x*x + c, t*x*y - s*z, t*x*z + s*y ],
    [ t*x*y + s*z, t*y*y + c, t*y*z - s*x ],
    [ t*x*z - s*y, t*y*z + s*x, t*z*z + c ]
  ]);
};













/**
 *  class FirminCSSMatrix
 *
 *  The [[FirminCSSMatrix]] class is a concrete implementation of the
 *  `CSSMatrix` interface defined in the [CSS 2D Transforms][2d] and
 *  [CSS 3D Transforms][3d] Module specifications.
 *
 *  [2d]: http://www.w3.org/TR/css3-2d-transforms/
 *  [3d]: http://www.w3.org/TR/css3-3d-transforms/
 *
 *  The implementation was largely copied from the `WebKitCSSMatrix` class, and
 *  the supparting maths libraries in the [WebKit][webkit] project. This is one
 *  reason why much of the code looks more like C++ than JavaScript.
 *
 *  [webkit]: http://webkit.org/
 *
 *  Its API is a superset of that provided by `WebKitCSSMatrix`, largely
 *  because various pieces of supporting code have been added as instance
 *  methods rather than pollute the global namespace. Examples of these include
 *  [[FirminCSSMatrix#isAffine]], [[FirminCSSMatrix#isIdentityOrTranslation]]
 *  and [[FirminCSSMatrix#adjoint]].
 **/

/**
 *  new FirminCSSMatrix(domstr)
 *  - domstr (String): a string representation of a 2D or 3D transform matrix
 *    in the form given by the CSS transform property, i.e. just like the
 *    output from [[FirminCSSMatrix#toString]].
 **/
FirminCSSMatrix = function(domstr) {
    this.m11 = this.m22 = this.m33 = this.m44 = 1;
    
               this.m12 = this.m13 = this.m14 =
    this.m21 =            this.m23 = this.m24 =
    this.m31 = this.m32 =            this.m34 =
    this.m41 = this.m42 = this.m43            = 0;
    
    if (typeof domstr == "string") {
        this.setMatrixValue(domstr);
    }
};

/**
 *  FirminCSSMatrix.displayName = "FirminCSSMatrix"
 **/
FirminCSSMatrix.displayName = "FirminCSSMatrix";

/**
 *  FirminCSSMatrix.degreesToRadians(angle) -> Number
 *  - angle (Number): an angle in degrees.
 *
 *  Converts angles in degrees, which are used by the external API, to angles
 *  in radians used in internal calculations.
 **/
FirminCSSMatrix.degreesToRadians = function(angle) {
    return angle * Math.PI / 180;
};

/**
 *  FirminCSSMatrix.determinant2x2(a, b, c, d) -> Number
 *  - a (Number): top-left value of the matrix.
 *  - b (Number): top-right value of the matrix.
 *  - c (Number): bottom-left value of the matrix.
 *  - d (Number): bottom-right value of the matrix.
 *
 *  Calculates the determinant of a 2x2 matrix.
 **/
FirminCSSMatrix.determinant2x2 = function(a, b, c, d) {
    return a * d - b * c;
};

/**
 *  FirminCSSMatrix.determinant3x3(matrix) -> Number
 *  - a1 (Number): matrix value in position [1, 1].
 *  - a2 (Number): matrix value in position [1, 2].
 *  - a3 (Number): matrix value in position [1, 3].
 *  - b1 (Number): matrix value in position [2, 1].
 *  - b2 (Number): matrix value in position [2, 2].
 *  - b3 (Number): matrix value in position [2, 3].
 *  - c1 (Number): matrix value in position [3, 1].
 *  - c2 (Number): matrix value in position [3, 2].
 *  - c3 (Number): matrix value in position [3, 3].
 *
 *  Calculates the determinant of a 3x3 matrix.
 **/
FirminCSSMatrix.determinant3x3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3) {
    var determinant2x2 = FirminCSSMatrix.determinant2x2;
    return a1 * determinant2x2(b2, b3, c2, c3) -
           b1 * determinant2x2(a2, a3, c2, c3) +
           c1 * determinant2x2(a2, a3, b2, b3);
};

/**
 *  FirminCSSMatrix.determinant4x4(matrix) -> Number
 *  - matrix (FirminCSSMatrix): the matrix to calculate the determinant of.
 *
 *  Calculates the determinant of a 4x4 matrix.
 **/
FirminCSSMatrix.determinant4x4 = function(m) {
    var determinant3x3 = FirminCSSMatrix.determinant3x3,
    
    // Assign to individual variable names to aid selecting correct elements
    a1 = m.m11, b1 = m.m21, c1 = m.m31, d1 = m.m41,
    a2 = m.m12, b2 = m.m22, c2 = m.m32, d2 = m.m42,
    a3 = m.m13, b3 = m.m23, c3 = m.m33, d3 = m.m43,
    a4 = m.m14, b4 = m.m24, c4 = m.m34, d4 = m.m44;
    
    return a1 * determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4) -
           b1 * determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4) +
           c1 * determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4) -
           d1 * determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
};

/**
 *  FirminCSSMatrix#a -> Number
 *  The first 2D vector value.
 **/

/**
 *  FirminCSSMatrix#b -> Number
 *  The second 2D vector value.
 **/

/**
 *  FirminCSSMatrix#c -> Number
 *  The third 2D vector value.
 **/

/**
 *  FirminCSSMatrix#d -> Number
 *  The fourth 2D vector value.
 **/

/**
 *  FirminCSSMatrix#e -> Number
 *  The fifth 2D vector value.
 **/

/**
 *  FirminCSSMatrix#f -> Number
 *  The sixth 2D vector value.
 **/

/**
 *  FirminCSSMatrix#m11 -> Number
 *  The 3D matrix value in the first row and first column.
 **/

/**
 *  FirminCSSMatrix#m12 -> Number
 *  The 3D matrix value in the first row and second column.
 **/

/**
 *  FirminCSSMatrix#m13 -> Number
 *  The 3D matrix value in the first row and third column.
 **/

/**
 *  FirminCSSMatrix#m14 -> Number
 *  The 3D matrix value in the first row and fourth column.
 **/

/**
 *  FirminCSSMatrix#m21 -> Number
 *  The 3D matrix value in the second row and first column.
 **/

/**
 *  FirminCSSMatrix#m22 -> Number
 *  The 3D matrix value in the second row and second column.
 **/

/**
 *  FirminCSSMatrix#m23 -> Number
 *  The 3D matrix value in the second row and third column.
 **/

/**
 *  FirminCSSMatrix#m24 -> Number
 *  The 3D matrix value in the second row and fourth column.
 **/

/**
 *  FirminCSSMatrix#m31 -> Number
 *  The 3D matrix value in the third row and first column.
 **/

/**
 *  FirminCSSMatrix#m32 -> Number
 *  The 3D matrix value in the third row and second column.
 **/

/**
 *  FirminCSSMatrix#m33 -> Number
 *  The 3D matrix value in the third row and third column.
 **/

/**
 *  FirminCSSMatrix#m34 -> Number
 *  The 3D matrix value in the third row and fourth column.
 **/

/**
 *  FirminCSSMatrix#m41 -> Number
 *  The 3D matrix value in the fourth row and first column.
 **/

/**
 *  FirminCSSMatrix#m42 -> Number
 *  The 3D matrix value in the fourth row and second column.
 **/

/**
 *  FirminCSSMatrix#m43 -> Number
 *  The 3D matrix value in the fourth row and third column.
 **/

/**
 *  FirminCSSMatrix#m44 -> Number
 *  The 3D matrix value in the fourth row and fourth column.
 **/

[["m11", "a"],
 ["m12", "b"],
 ["m21", "c"],
 ["m22", "d"],
 ["m41", "e"],
 ["m42", "f"]].forEach(function(pair) {
    var key3d = pair[0], key2d = pair[1];
    
    Object.defineProperty(FirminCSSMatrix.prototype, key2d, {
        set: function(val) {
            this[key3d] = val;
        },
        
        get: function() {
            return this[key3d];
        }
    });
});

/**
 *  FirminCSSMatrix#isAffine() -> Boolean
 *
 *  Determines whether the matrix is affine.
 **/
FirminCSSMatrix.prototype.isAffine = function() {
    return this.m13 === 0 && this.m14 === 0 &&
           this.m23 === 0 && this.m24 === 0 &&
           this.m31 === 0 && this.m32 === 0 &&
           this.m33 === 1 && this.m34 === 0 &&
           this.m43 === 0 && this.m44 === 1;
};

/**
 *  FirminCSSMatrix#multiply(otherMatrix) -> FirminCSSMatrix
 *  - otherMatrix (FirminCSSMatrix): the matrix to multiply this one by.
 *
 *  Multiplies the matrix by a given matrix and returns the result.
 **/
FirminCSSMatrix.prototype.multiply = function(otherMatrix) {
    var a = otherMatrix,
        b = this,
        c = new FirminCSSMatrix();
    
    c.m11 = a.m11 * b.m11 + a.m12 * b.m21 + a.m13 * b.m31 + a.m14 * b.m41;
    c.m12 = a.m11 * b.m12 + a.m12 * b.m22 + a.m13 * b.m32 + a.m14 * b.m42;
    c.m13 = a.m11 * b.m13 + a.m12 * b.m23 + a.m13 * b.m33 + a.m14 * b.m43;
    c.m14 = a.m11 * b.m14 + a.m12 * b.m24 + a.m13 * b.m34 + a.m14 * b.m44;
    
    c.m21 = a.m21 * b.m11 + a.m22 * b.m21 + a.m23 * b.m31 + a.m24 * b.m41;
    c.m22 = a.m21 * b.m12 + a.m22 * b.m22 + a.m23 * b.m32 + a.m24 * b.m42;
    c.m23 = a.m21 * b.m13 + a.m22 * b.m23 + a.m23 * b.m33 + a.m24 * b.m43;
    c.m24 = a.m21 * b.m14 + a.m22 * b.m24 + a.m23 * b.m34 + a.m24 * b.m44;
    
    c.m31 = a.m31 * b.m11 + a.m32 * b.m21 + a.m33 * b.m31 + a.m34 * b.m41;
    c.m32 = a.m31 * b.m12 + a.m32 * b.m22 + a.m33 * b.m32 + a.m34 * b.m42;
    c.m33 = a.m31 * b.m13 + a.m32 * b.m23 + a.m33 * b.m33 + a.m34 * b.m43;
    c.m34 = a.m31 * b.m14 + a.m32 * b.m24 + a.m33 * b.m34 + a.m34 * b.m44;
    
    c.m41 = a.m41 * b.m11 + a.m42 * b.m21 + a.m43 * b.m31 + a.m44 * b.m41;
    c.m42 = a.m41 * b.m12 + a.m42 * b.m22 + a.m43 * b.m32 + a.m44 * b.m42;
    c.m43 = a.m41 * b.m13 + a.m42 * b.m23 + a.m43 * b.m33 + a.m44 * b.m43;
    c.m44 = a.m41 * b.m14 + a.m42 * b.m24 + a.m43 * b.m34 + a.m44 * b.m44;
    
    return c;
};

/**
 *  FirminCSSMatrix#isIdentityOrTranslation() -> Boolean
 *
 *  Returns whether the matrix is the identity matrix or a translation matrix.
 **/
FirminCSSMatrix.prototype.isIdentityOrTranslation = function() {
    var t = this;
    return t.m11 === 1 && t.m12 === 0 && t.m13 === 0 && t.m14 === 0 &&
           t.m21 === 0 && t.m22 === 1 && t.m23 === 0 && t.m24 === 0 &&
           t.m31 === 0 && t.m31 === 0 && t.m33 === 1 && t.m34 === 0 &&
    /* m41, m42 and m43 are the translation points */   t.m44 === 1;
};

/**
 *  FirminCSSMatrix#adjoint() -> FirminCSSMatrix
 *
 *  Returns the adjoint matrix.
 **/
FirminCSSMatrix.prototype.adjoint = function() {
    var result = new FirminCSSMatrix(), t = this,
        determinant3x3 = FirminCSSMatrix.determinant3x3,
        
        a1 = t.m11, b1 = t.m12, c1 = t.m13, d1 = t.m14,
        a2 = t.m21, b2 = t.m22, c2 = t.m23, d2 = t.m24,
        a3 = t.m31, b3 = t.m32, c3 = t.m33, d3 = t.m34,
        a4 = t.m41, b4 = t.m42, c4 = t.m43, d4 = t.m44;
    
    // Row column labeling reversed since we transpose rows & columns
    result.m11 =  determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
    result.m21 = -determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
    result.m31 =  determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
    result.m41 = -determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
    
    result.m12 = -determinant3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
    result.m22 =  determinant3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
    result.m32 = -determinant3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
    result.m42 =  determinant3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);
    
    result.m13 =  determinant3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
    result.m23 = -determinant3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
    result.m33 =  determinant3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
    result.m43 = -determinant3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);
    
    result.m14 = -determinant3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
    result.m24 =  determinant3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
    result.m34 = -determinant3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
    result.m44 =  determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
    
    return result;
};

/**
 *  FirminCSSMatrix#inverse() -> FirminCSSMatrix | null
 *
 *  If the matrix is invertible, returns its inverse, otherwise returns null.
 **/
FirminCSSMatrix.prototype.inverse = function() {
    var inv, det, result, i, j;
    
    if (this.isIdentityOrTranslation()) {
        inv = new FirminCSSMatrix();
        
        if (!(this.m41 === 0 && this.m42 === 0 && this.m43 === 0)) {
            inv.m41 = -this.m41;
            inv.m42 = -this.m42;
            inv.m43 = -this.m43;
        }
        
        return inv;
    }
    
    // Calculate the adjoint matrix
    result = this.adjoint();
    
    // Calculate the 4x4 determinant
    det = FirminCSSMatrix.determinant4x4(this);
    
    // If the determinant is zero, then the inverse matrix is not unique
    if (Math.abs(det) < 1e-8) return null;
    
    // Scale the adjoint matrix to get the inverse
    for (i = 1; i < 5; i++) {
        for (j = 1; j < 5; j++) {
            result[("m" + i) + j] /= det;
        }
    }
    
    return result;
};

/**
 *  FirminCSSMatrix#rotate(rotX, rotY, rotZ) -> FirminCSSMatrix
 *  - rotX (Number): the rotation around the x axis.
 *  - rotY (Number): the rotation around the y axis. If undefined, the x
 *    component is used.
 *  - rotZ (Number): the rotation around the z axis. If undefined, the x
 *    component is used.
 *
 *  Returns the result of rotating the matrix by a given vector.
 *
 *  If only the first argument is provided, the matrix is only rotated about
 *  the z axis.
 **/
FirminCSSMatrix.prototype.rotate = function(rx, ry, rz) {
    var degreesToRadians = FirminCSSMatrix.degreesToRadians;
    
    if (typeof rx != "number" || isNaN(rx)) rx = 0;
    
    if ((typeof ry != "number" || isNaN(ry)) &&
        (typeof rz != "number" || isNaN(rz))) {
        rz = rx;
        rx = 0;
        ry = 0;
    }
    
    if (typeof ry != "number" || isNaN(ry)) ry = 0;
    if (typeof rz != "number" || isNaN(rz)) rz = 0;
    
    rx = degreesToRadians(rx);
    ry = degreesToRadians(ry);
    rz = degreesToRadians(rz);
    
    var tx = new FirminCSSMatrix(),
        ty = new FirminCSSMatrix(),
        tz = new FirminCSSMatrix(),
        sinA, cosA, sinA2;
    
    rz /= 2;
    sinA = Math.sin(rz);
    cosA = Math.cos(rz);
    sinA2 = sinA * sinA;
    
    // Matrices are identity outside the assigned values
    tz.m11 = tz.m22 = 1 - 2 * sinA2;
    tz.m12 = tz.m21 = 2 * sinA * cosA;
    tz.m21 *= -1;
    
    ry /= 2;
    sinA  = Math.sin(ry);
    cosA  = Math.cos(ry);
    sinA2 = sinA * sinA;
    
    ty.m11 = ty.m33 = 1 - 2 * sinA2;
    ty.m13 = ty.m31 = 2 * sinA * cosA;
    ty.m13 *= -1;
    
    rx /= 2;
    sinA = Math.sin(rx);
    cosA = Math.cos(rx);
    sinA2 = sinA * sinA;
    
    tx.m22 = tx.m33 = 1 - 2 * sinA2;
    tx.m23 = tx.m32 = 2 * sinA * cosA;
    tx.m32 *= -1;
    
    return tz.multiply(ty).multiply(tx).multiply(this);
};

/**
 *  FirminCSSMatrix#rotateAxisAngle(rotX, rotY, rotZ, angle) -> FirminCSSMatrix
 *  - rotX (Number): the rotation around the x axis.
 *  - rotY (Number): the rotation around the y axis. If undefined, the x
 *    component is used.
 *  - rotZ (Number): the rotation around the z axis. If undefined, the x
 *    component is used.
 *  - angle (Number): the angle of rotation about the axis vector, in degrees.
 *
 *  Returns the result of rotating the matrix around a given vector by a given
 *  angle.
 *
 *  If the given vector is the origin vector then the matrix is rotated by the
 *  given angle around the z axis.
 **/
FirminCSSMatrix.prototype.rotateAxisAngle = function(x, y, z, a) {
    if (typeof x != "number" || isNaN(x)) x = 0;
    if (typeof y != "number" || isNaN(y)) y = 0;
    if (typeof z != "number" || isNaN(z)) z = 0;
    if (typeof a != "number" || isNaN(a)) a = 0;
    if (x === 0 && y === 0 && z === 0) z = 1;
    
    var t   = new FirminCSSMatrix(),
        len = Math.sqrt(x * x + y * y + z * z),
        cosA, sinA, sinA2, csA, x2, y2, z2;
    
    a     = (FirminCSSMatrix.degreesToRadians(a) || 0) / 2;
    cosA  = Math.cos(a);
    sinA  = Math.sin(a);
    sinA2 = sinA * sinA;
    
    // Bad vector, use something sensible
    if (len === 0) {
        x = 0;
        y = 0;
        z = 1;
    } else if (len !== 1) {
        x /= len;
        y /= len;
        z /= len;
    }
    
    // Optimise cases where axis is along major axis
    if (x === 1 && y === 0 && z === 0) {
        t.m22 = t.m33 = 1 - 2 * sinA2;
        t.m23 = t.m32 = 2 * cosA * sinA;
        t.m32 *= -1;
    } else if (x === 0 && y === 1 && z === 0) {
        t.m11 = t.m33 = 1 - 2 * sinA2;
        t.m13 = t.m31 = 2 * cosA * sinA;
        t.m13 *= -1;
    } else if (x === 0 && y === 0 && z === 1) {
        t.m11 = t.m22 = 1 - 2 * sinA2;
        t.m12 = t.m21 = 2 * cosA * sinA;
        t.m21 *= -1;
    } else {
        csA = sinA * cosA;
        x2  = x * x;
        y2  = y * y;
        z2  = z * z;
        
        t.m11 = 1 - 2 * (y2 + z2) * sinA2;
        t.m12 = 2 * (x * y * sinA2 + z * csA);
        t.m13 = 2 * (x * z * sinA2 - y * csA);
        t.m21 = 2 * (y * x * sinA2 - z * csA);
        t.m22 = 1 - 2 * (z2 + x2) * sinA2;
        t.m23 = 2 * (y * z * sinA2 + x * csA);
        t.m31 = 2 * (z * x * sinA2 + y * csA);
        t.m32 = 2 * (z * y * sinA2 - x * csA);
        t.m33 = 1 - 2 * (x2 + y2) * sinA2;
    }
    
    return this.multiply(t);
};

/**
 *  FirminCSSMatrix#scale(scaleX, scaleY, scaleZ) -> FirminCSSMatrix
 *  - scaleX (Number): the scaling factor in the x axis.
 *  - scaleY (Number): the scaling factor in the y axis. If undefined, the x
 *    component is used.
 *  - scaleZ (Number): the scaling factor in the z axis. If undefined, 1 is
 *    used.
 *
 *  Returns the result of scaling the matrix by a given vector.
 **/
FirminCSSMatrix.prototype.scale = function(scaleX, scaleY, scaleZ) {
    var transform = new FirminCSSMatrix();
    
    if (typeof scaleX != "number" || isNaN(scaleX)) scaleX = 1;
    if (typeof scaleY != "number" || isNaN(scaleY)) scaleY = scaleX;
    if (typeof scaleZ != "number" || isNaN(scaleZ)) scaleZ = 1;
    
    transform.m11 = scaleX;
    transform.m22 = scaleY;
    transform.m33 = scaleZ;
    
    return this.multiply(transform);
};

/**
 *  FirminCSSMatrix#translate(x, y, z) -> FirminCSSMatrix
 *  - x (Number): the x component of the vector.
 *  - y (Number): the y component of the vector.
 *  - z (Number): the z component of the vector. If undefined, 0 is used.
 *
 *  Returns the result of translating the matrix by a given vector.
 **/
FirminCSSMatrix.prototype.translate = function(x, y, z) {
    var t = new FirminCSSMatrix();
    
    if (typeof x != "number" || isNaN(x)) x = 0;
    if (typeof y != "number" || isNaN(y)) y = 0;
    if (typeof z != "number" || isNaN(z)) z = 0;
    
    t.m41 = x;
    t.m42 = y;
    t.m43 = z;
    
    return this.multiply(t);
};

/**
 *  FirminCSSMatrix#setMatrixValue(domstr) -> undefined
 *  - domstr (String): a string representation of a 2D or 3D transform matrix
 *    in the form given by the CSS transform property, i.e. just like the
 *    output from [[FirminCSSMatrix#toString]].
 *
 *  Sets the matrix values using a string representation, such as that produced
 *  by the [[FirminCSSMatrix#toString]] method.
 **/
FirminCSSMatrix.prototype.setMatrixValue = function(domstr) {
        domstr = domstr.trim();
    var mstr   = domstr.match(/^matrix(3d)?\(\s*(.+)\s*\)$/),
        is3d, chunks, len, points, i, chunk;
    
    if (!mstr) return;
    
    is3d   = !!mstr[1];
    chunks = mstr[2].split(/\s*,\s*/);
    len    = chunks.length;
    points = new Array(len);
    
    if ((is3d && len !== 16) || !(is3d || len === 6)) return;
    
    for (i = 0; i < len; i++) {
        chunk = chunks[i];
        if (chunk.match(/^-?\d+(\.\d+)?$/)) {
            points[i] = parseFloat(chunk);
        } else return;
    }
    
    for (i = 0; i < len; i++) {
        point = is3d ?
            ("m" + (Math.floor(i / 4) + 1)) + (i % 4 + 1) :
            String.fromCharCode(i + 97); // ASCII char 97 == 'a'
        this[point] = points[i];
    }
};

/**
 *  FirminCSSMatrix#toString() -> String
 *
 *  Returns a string representation of the matrix.
 **/
FirminCSSMatrix.prototype.toString = function() {
    var self = this, points, prefix;
    
    if (this.isAffine()) {
        prefix = "matrix(";
        points = ["a", "b", "c", "d", "e", "f"];
    } else {
        prefix = "matrix3d(";
        points = ["m11", "m12", "m13", "m14",
                  "m21", "m22", "m23", "m24",
                  "m31", "m32", "m33", "m34",
                  "m41", "m42", "m43", "m44"];
    }
    
    return prefix + points.map(function(p) {
        return self[p].toFixed(6);
    }).join(", ") + ")";
};














/*
 * @preserve Morf v0.1.5
 * http://www.joelambert.co.uk/morf
 *
 * Copyright 2011, Joe Lambert.
 * Free to use under the MIT license.
 * http://www.opensource.org/licenses/mit-license.php
 */

var CSSMatrixDecomposed = function(obj) {
	obj === undefined ? obj = {} : null;
	var components = {perspective: null, translate: null, skew: null, scale: null, rotate: null};
	
	for(var i in components)
		this[i] = obj[i] ? obj[i] : new Vector4();

	/**
	 * Tween between two decomposed matrices
	 * @param {CSSMatrixDecomposed} dm The destination decomposed matrix
	 * @param {float} progress A float value between 0-1, representing the percentage of completion
	 * @param {function} fn An easing function following the prototype function(pos){}
	 * @author Joe Lambert
	 * @returns {WebKitCSSMatrix} A new matrix for the tweened state
	 */
		
	this.tween = function(dm, progress, fn) {
		if(fn === undefined)
			fn = function(pos) {return pos;}; // Default to a linear easing
		
		if(!dm)
			dm = new CSSMatrixDecomposed(new FirminCSSMatrix().decompose());
		
		var r = new CSSMatrixDecomposed(),
			i = index = null,
			trans = '';
		
		progress = fn(progress);

		for(index in components)
			for(i in {x:'x', y:'y', z:'z', w:'w'})
				r[index][i] = (this[index][i] + (dm[index][i] - this[index][i]) * progress ).toFixed(5);

		trans = 'matrix3d(1,0,0,0, 0,1,0,0, 0,0,1,0, '+r.perspective.x+', '+r.perspective.y+', '+r.perspective.z+', '+r.perspective.w+') ' +
				'translate3d('+r.translate.x+'px, '+r.translate.y+'px, '+r.translate.y+'px) ' +
				'rotateX('+r.rotate.x+'rad) rotateY('+r.rotate.y+'rad) rotateZ('+r.rotate.z+'rad) ' +
				'matrix3d(1,0,0,0, 0,1,0,0, 0,'+r.skew.z+',1,0, 0,0,0,1) ' +
				'matrix3d(1,0,0,0, 0,1,0,0, '+r.skew.y+',0,1,0, 0,0,0,1) ' +
				'matrix3d(1,0,0,0, '+r.skew.x+',1,0,0, 0,0,1,0, 0,0,0,1) ' +
				'scale3d('+r.scale.x+', '+r.scale.y+', '+r.scale.z+')';

		try { r = new FirminCSSMatrix(trans); return r; }
		catch(e) { console.error('Invalid matrix string: '+trans); return '' };
	};
};

var Vector4 = function(x, y, z, w)
{
	this.x = x ? x : 0;
	this.y = y ? y : 0;
	this.z = z ? z : 0;
	this.w = w ? w : 0;
	
	
	/**
	 * Ensure that values are not undefined
	 * @author Joe Lambert
	 * @returns null
	 */
	
	this.checkValues = function() {
		this.x = this.x ? this.x : 0;
		this.y = this.y ? this.y : 0;
		this.z = this.z ? this.z : 0;
		this.w = this.w ? this.w : 0;
	};
	
	
	/**
	 * Get the length of the vector
	 * @author Joe Lambert
	 * @returns {float}
	 */
	
	this.length = function() {
		this.checkValues();
		return Math.sqrt(this.x*this.x + this.y*this.y + this.z*this.z);
	};
	
	
	/**
	 * Get a normalised representation of the vector
	 * @author Joe Lambert
	 * @returns {Vector4}
	 */
	
	this.normalise = function() {
		var len = this.length(),
			v = new Vector4(this.x / len, this.y / len, this.z / len);
		
		return v;
	};


	/**
	 * Vector Dot-Product
	 * @param {Vector4} v The second vector to apply the product to
	 * @author Joe Lambert
	 * @returns {float} The Dot-Product of this and v.
	 */

	this.dot = function(v) {
		return this.x*v.x + this.y*v.y + this.z*v.z + this.w*v.w;
	};
	
	
	/**
	 * Vector Cross-Product
	 * @param {Vector4} v The second vector to apply the product to
	 * @author Joe Lambert
	 * @returns {Vector4} The Cross-Product of this and v.
	 */
	
	this.cross = function(v) {
		return new Vector4(this.y*v.z - this.z*v.y, this.z*v.x - this.x*v.z, this.x*v.y - this.y*v.x);
	};
	

	/**
	 * Helper function required for matrix decomposition
	 * A Javascript implementation of pseudo code available from http://www.w3.org/TR/css3-2d-transforms/#matrix-decomposition
	 * @param {Vector4} aPoint A 3D point
	 * @param {float} ascl 
	 * @param {float} bscl
	 * @author Joe Lambert
	 * @returns {Vector4}
	 */
	
	this.combine = function(aPoint, ascl, bscl) {
		return new Vector4( (ascl * this.x) + (bscl * aPoint.x), 
							(ascl * this.y) + (bscl * aPoint.y), 
							(ascl * this.z) + (bscl * aPoint.z) );
	}
};

FirminCSSMatrix.prototype.determinant = function() {
	return 	this.m14 * this.m23 * this.m32 * this.m41-this.m13 * this.m24 * this.m32 * this.m41 -
			this.m14 * this.m22 * this.m33 * this.m41+this.m12 * this.m24 * this.m33 * this.m41 +
			this.m13 * this.m22 * this.m34 * this.m41-this.m12 * this.m23 * this.m34 * this.m41 -
			this.m14 * this.m23 * this.m31 * this.m42+this.m13 * this.m24 * this.m31 * this.m42 +
			this.m14 * this.m21 * this.m33 * this.m42-this.m11 * this.m24 * this.m33 * this.m42 -
			this.m13 * this.m21 * this.m34 * this.m42+this.m11 * this.m23 * this.m34 * this.m42 +
			this.m14 * this.m22 * this.m31 * this.m43-this.m12 * this.m24 * this.m31 * this.m43 -
			this.m14 * this.m21 * this.m32 * this.m43+this.m11 * this.m24 * this.m32 * this.m43 +
			this.m12 * this.m21 * this.m34 * this.m43-this.m11 * this.m22 * this.m34 * this.m43 -
			this.m13 * this.m22 * this.m31 * this.m44+this.m12 * this.m23 * this.m31 * this.m44 +
			this.m13 * this.m21 * this.m32 * this.m44-this.m11 * this.m23 * this.m32 * this.m44 -
			this.m12 * this.m21 * this.m33 * this.m44+this.m11 * this.m22 * this.m33 * this.m44;
};

FirminCSSMatrix.prototype.decompose = function() {
	var matrix = new FirminCSSMatrix(this.toString()),
		perspectiveMatrix = rightHandSide = inversePerspectiveMatrix = transposedInversePerspectiveMatrix =
		perspective = translate = row = i = scale = skew = pdum3 =  rotate = null;
	
	if (matrix.m33 == 0)
	    return new CSSMatrixDecomposed(new FirminCSSMatrix().decompose()); // Return the identity matrix

	// Normalize the matrix.
	for (i = 1; i <= 4; i++)
	    for (j = 1; j <= 4; j++)
	        matrix['m'+i+j] /= matrix.m44;

	// perspectiveMatrix is used to solve for perspective, but it also provides
	// an easy way to test for singularity of the upper 3x3 component.
	perspectiveMatrix = matrix;

	for (i = 1; i <= 3; i++)
	    perspectiveMatrix['m'+i+'4'] = 0;

	perspectiveMatrix.m44 = 1;

	if (perspectiveMatrix.determinant() == 0)
	    return new CSSMatrixDecomposed(new FirminCSSMatrix().decompose()); // Return the identity matrix

	// First, isolate perspective.
	if (matrix.m14 != 0 || matrix.m24 != 0 || matrix.m34 != 0)
	{
	    // rightHandSide is the right hand side of the equation.
		rightHandSide = new Vector4(matrix.m14, matrix.m24, matrix.m34, matrix.m44);
		
	    // Solve the equation by inverting perspectiveMatrix and multiplying
	    // rightHandSide by the inverse.
	    inversePerspectiveMatrix 			= perspectiveMatrix.inverse();
	    transposedInversePerspectiveMatrix 	= inversePerspectiveMatrix.transpose();
	    perspective 						= transposedInversePerspectiveMatrix.transformVector(rightHandSide);

	     // Clear the perspective partition
	    matrix.m14 = matrix.m24 = matrix.m34 = 0;
	    matrix.m44 = 1;
	}
	else
	{
		// No perspective.
		perspective = new Vector4(0,0,0,1);
	}

	// Next take care of translation
	translate = new Vector4(matrix.m41, matrix.m42, matrix.m43);

	matrix.m41 = 0;
	matrix.m42 = 0;
	matrix.m43 = 0;	
	
	// Now get scale and shear. 'row' is a 3 element array of 3 component vectors
	row = [
		new Vector4(), new Vector4(), new Vector4()
	];
	
	for (i = 1; i <= 3; i++)
	{
		row[i-1].x = matrix['m'+i+'1'];
	    row[i-1].y = matrix['m'+i+'2'];
	    row[i-1].z = matrix['m'+i+'3'];
	}

	// Compute X scale factor and normalize first row.
	scale = new Vector4();
	skew = new Vector4();
	
	scale.x = row[0].length();
	row[0] = row[0].normalise();
	
	// Compute XY shear factor and make 2nd row orthogonal to 1st.
	skew.x = row[0].dot(row[1]);
	row[1] = row[1].combine(row[0], 1.0, -skew.x);
	
	// Now, compute Y scale and normalize 2nd row.
	scale.y = row[1].length();
	row[1] = row[1].normalise();
	skew.x /= scale.y;
	
	// Compute XZ and YZ shears, orthogonalize 3rd row
	skew.y = row[0].dot(row[2]);
	row[2] = row[2].combine(row[0], 1.0, -skew.y);
	skew.z = row[1].dot(row[2]);
	row[2] = row[2].combine(row[1], 1.0, -skew.z);
	
	// Next, get Z scale and normalize 3rd row.
	scale.z = row[2].length();
	row[2] = row[2].normalise();
	skew.y /= scale.z;
	skew.y /= scale.z;
	
	// At this point, the matrix (in rows) is orthonormal.
	// Check for a coordinate system flip.  If the determinant
	// is -1, then negate the matrix and the scaling factors.
	pdum3 = row[1].cross(row[2])
	if (row[0].dot(pdum3) < 0)
	{
		for (i = 0; i < 3; i++)
		{
	        scale.x *= -1;
	        row[i].x *= -1;
	        row[i].y *= -1;
	        row[i].z *= -1;	
		}
	}

	// Now, get the rotations out
	rotate = new Vector4();
	rotate.y = Math.asin(-row[0].z);
	if (Math.cos(rotate.y) != 0)
	{
		rotate.x = Math.atan2(row[1].z, row[2].z);
		rotate.z = Math.atan2(row[0].y, row[0].x);
	}
	else
	{
		rotate.x = Math.atan2(-row[2].x, row[1].y);
		rotate.z = 0;
	}
	
	return new CSSMatrixDecomposed({
		perspective: perspective,
		translate: translate,
		skew: skew,
		scale: scale,
		rotate: rotate
	});
};