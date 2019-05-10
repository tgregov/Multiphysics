// 2D Rectangle
// List of points
d = 0.4;
x = 10;
y = 0.1;
Point(1) = {0, 0, 0, d};
Point(2) = {x, 0, 0, d};
Point(3) = {0, y, 0, d};
Point(4) = {x, y, 0, d};

// List of lines
Line(1) = {3, 1};
Line(2) = {4, 3};
Line(3) = {2, 4};
Line(4) = {1, 2};

// Surface
Curve Loop(1) = {1, 4, 3, 2};
Plane Surface(1) = {1};

Physical Surface(1) = {1};
// Boundary conditions
Physical Curve("BC_Left") = {1};
Physical Curve("BC_Up") = {2};
Physical Curve("BC_Down") = {4};
Physical Curve("BC_Right") = {3};

Transfinite Surface{1};