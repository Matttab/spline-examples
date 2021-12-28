# spline-examples
Use of simple polynominal and B-Splines in GNU Octave and Matlab for geometrical shape approximation

The objective are simple spline approximations of common geometric shapes with potential practical use in CNC machining.
Interest is on spline interpolations with a low number of parameters, describing common shapes such as circles, ellipses,
squares, slots, pockets and contours of arbitrary shape.

Splines are used in the CAD to define the geometry, in the CAM to define the tool path by
a list of dedicated G-code machine instructions, and in the CNC control SW, to eventually
interpolate G5 spline instruction to actual workpiece or tool movements.

Mapping of spline parameters to dimensions of the geometrical entities is welcome.
Conversion of parameters between different types of splines is part of adaption to the 
broad varyity of commercial and home-made hobbist CAD-CAM-CNC systems.
