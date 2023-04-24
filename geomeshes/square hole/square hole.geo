SetFactory("OpenCASCADE");

// WALLS OF ROOM
top = 10;
bottom = 0;
left = 0;
right = 10;
// corners
Point(1) = {right, bottom, 0};
Point(2) = {right, top, 0};
Point(3) = {left, top, 0};
Point(4) = {left, bottom, 0};
// edges
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// loop
Curve Loop(1) = {1, 2, 3, 4};
Physical Curve("Walls", 1) = {1, 2, 3, 4};

// SPEAKER BOX
spkr_t = 1.4;
spkr_b = 0.6;
spkr_bin = 0.7;
spkr_l = 4.6;
spkr_lin = 4.7;
spkr_r = 5.4;
spkr_rin = 5.3;
// corners
Point(5) = {spkr_r, spkr_b, 0};
Point(6) = {spkr_r, spkr_t, 0};
Point(7) = {spkr_rin, spkr_t, 0};
Point(8) = {spkr_rin, spkr_bin, 0};
Point(9) = {spkr_lin, spkr_bin, 0};
Point(10) = {spkr_lin, spkr_t, 0};
Point(11) = {spkr_l, spkr_t, 0};
Point(12) = {spkr_l, spkr_b, 0};
// edges
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 5};
// loop
Curve Loop(3) = {6, 7, 8, 9, 10, 11, 12, 13}; // speaker
Physical Curve("Speaker", 2) = {6, 7, 8, 9, 10, 11, 12, 13};

// SOURCE CIRCLE
Circle(5) = {5, 1, 0, 0.2, 0, 2*Pi};
Curve Loop(2) = {5}; // source

// HOLE
hole_t = 7.5;
hole_b = 2.5;
hole_l = 0.5;
hole_r = 9.5;
// corners
Point(14) = {hole_r, hole_b, 0};
Point(15) = {hole_r, hole_t, 0};
Point(16) = {hole_l, hole_t, 0};
Point(17) = {hole_l, hole_b, 0};
// edges
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {17, 14};
// loop
Curve Loop(4) = {14, 15, 16, 17};
Physical Curve("Hole", 3) = {4};

// SURFACES
Plane Surface(1) = {1, 2, 3, 4}; // room surface
Plane Surface(2) = {2}; // source surface

Physical Surface("Room", 1) = {1};
Physical Surface("Source", 2) = {2};
