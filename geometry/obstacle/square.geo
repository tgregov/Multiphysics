// Square obstacle
// List of points
d = 0.03;
x = 1;
y = 1;
xC = 0.45;
yC = 0.45;
LxC = 0.1;
LyC = 0.1;

Point(1) = {0, 0, 0, d};
Point(2) = {x, 0, 0, d};
Point(3) = {0, y, 0, d};
Point(4) = {x, y, 0, d};
Point(5) = {xC, yC, 0, d};
Point(6) = {xC + LxC, yC, 0, d};
Point(7) = {xC, yC + LyC, 0, d};
Point(8) = {xC + LxC, yC + LyC, 0, d};

// List of lines
Line(1) = {3, 1};
Line(2) = {4, 3};
Line(3) = {2, 4};
Line(4) = {1, 2};
Line(5) = {7, 5};
Line(6) = {8, 7};
Line(7) = {6, 8};
Line(8) = {5, 6};

// Surface
Curve Loop(1) = {1, 4, 3, 2};
Curve Loop(2) = {5, 8, 7, 6};
Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};

// Boundary conditions
Physical Curve("BC_Left") = {1};
Physical Curve("BC_Up") = {2};
Physical Curve("BC_Down") = {4};
Physical Curve("BC_Right") = {3};

Physical Curve("BC_Obstacle") = {5, 6, 7 ,8};

// Transfinite Surface{1};

