//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 10, 0, 1.0};
//+
Point(4) = {0, 10, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Outer Walls", 1) = {1, 2, 3, 4};
//+
Physical Surface("Room", 1) = {1};
