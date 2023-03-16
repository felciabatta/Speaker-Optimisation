// Gmsh project created on Thu Mar 16 16:16:16 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Physical Curve("end wall") = {1};
//+
Physical Curve("front wall") = {3};
//+
Physical Curve("side wall") = {4, 2};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Surface("The surface") = {1};
//+
Physical Surface("The surface") += {1};
//+
Physical Surface("The surface") += {1};
//+
Physical Surface("The surface") += {1};
