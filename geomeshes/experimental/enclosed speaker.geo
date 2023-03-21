SetFactory("OpenCASCADE");

// corners of room
top = 10;
bottom = 0;
left = 0;
right = 10;
Point(1) = {right, bottom, 0};
Point(2) = {right, top, 0};
Point(3) = {left, top, 0};
Point(4) = {left, bottom, 0};

// corners of speaker box
spkr_t = 1.4;
spkr_tin = 1.3;
spkr_b = 0.6;
spkr_bin = 0.7;
spkr_l = 4.6;
spkr_lin = 4.7;
spkr_r = 5.4;
spkr_rin = 5.3;
// speaker out
Point(5) = {spkr_r, spkr_b, 0};
Point(6) = {spkr_r, spkr_t, 0};
Point(7) = {spkr_l, spkr_t, 0};
Point(8) = {spkr_l, spkr_b, 0};
// speaker in
Point(9) = {spkr_rin, spkr_bin, 0};
Point(10) = {spkr_lin, spkr_bin, 0};
Point(11) = {spkr_lin, spkr_tin, 0};
Point(12) = {spkr_rin, spkr_tin, 0};

// walls of room
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Physical Curve("Walls", 1) = {1, 2, 3, 4};

// source circle
Circle(5) = {5, 1, 0, 0.2, 0, 2*Pi};

// speaker sides
// outer
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};
Line(9) = {8, 5};
// inner
Line(10) = {9, 10};
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 9};
Physical Curve("Speaker", 2) = {6, 7, 8, 9, 10, 11, 12, 13};

// define surface boundaries
Curve Loop(1) = {1, 2, 3, 4}; // walls
Curve Loop(2) = {5}; // source
Curve Loop(3) = {6, 7, 8, 9}; // speaker out
Curve Loop(4) = {10, 11, 12, 13}; // speaker in

// define surfaces
Plane Surface(1) = {1, 3}; // room surface
Physical Surface("Room", 1) = {1};

Plane Surface(2) = {2}; // source surface
Physical Surface("Source", 2) = {2};

Plane Surface(3) = {4, 2}; // inside speaker
Physical Surface("In Speaker", 3) = {3};
