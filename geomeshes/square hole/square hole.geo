//+
Point(1) = {0, -0, 0, 1.0};
//+
Point(2) = {0, 10, 0, 1.0};
//+
Point(3) = {10, 10, 0, 1.0};
//+
Point(4) = {10, 0, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Point(5) = {2.5, 7.5, 0, 1.0};
//+
Point(6) = {2.5, 2.5, 0, 1.0};
//+
Point(7) = {7.5, 7.5, 0, 1.0};
//+
Point(8) = {7.5, 2.5, 0, 1.0};
//+
Line(5) = {5, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 6};
//+
Line(8) = {6, 5};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Curve("Outer Walls", 1) = {1, 4, 3, 2};
//+
Physical Curve("Inner Walls", 2) = {5, 8, 7, 6};
//+
Physical Surface("Outer Room", 1) = {1};
//+
Physical Surface("Inner Room", 2) = {2};
