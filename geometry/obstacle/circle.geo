// Square obstacle
// List of points
d = 0.03;
dBis = 0.5;
x = 1;
y = 1;
xC = 0.5;
yC = 0.5;

Point(1) = {0, 0, 0, d};
Point(2) = {x, 0, 0, d};
Point(3) = {0, y, 0, d};
Point(4) = {x, y, 0, d};

Point(5) = {xC+0.075,	yC+0, 0, dBis};
Point(6) = {xC+0.0693, 	yC+0.0287, 0, dBis};
Point(7) = {xC+0.0530, 	yC+0.053, 0, dBis};
Point(8) = {xC+0.0287, 	yC+0.0693, 0, dBis};
Point(9) = {xC+0, 	yC+0.075, 0, dBis};
Point(10) = {xC-0.0287, yC+0.0693, 0, dBis};
Point(11) = {xC-0.0530, yC+0.053, 0, dBis};
Point(12) = {xC-0.0693, yC+0.0287, 0, dBis};
Point(13) = {xC-0.075, 	yC+0, 0, dBis};
Point(14) = {xC-0.0693, yC-0.0287, 0, dBis};
Point(15) = {xC-0.0530, yC-0.053, 0, dBis};
Point(16) = {xC-0.0287, yC-0.0693, 0, dBis};
Point(17) = {xC-0, 	yC-0.075, 0, dBis};
Point(18) = {xC+0.0287, yC-0.0693, 0, dBis};
Point(19) = {xC+0.0530, yC-0.053, 0, dBis};
Point(20) = {xC+0.0693, yC-0.0287, 0, dBis};

// List of lines
Line(1) = {3, 1};
Line(2) = {4, 3};
Line(3) = {2, 4};
Line(4) = {1, 2};
Line(5) = {5:20,5};

// Surface
Curve Loop(1) = {1, 4, 3, 2};
Curve Loop(2) = {5};
Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};

// Boundary conditions
Physical Curve("BC_Left") = {1};
Physical Curve("BC_Up") = {2};
Physical Curve("BC_Down") = {4};
Physical Curve("BC_Right") = {3};

Physical Curve("BC_Obstacle") = {5};

// Transfinite Surface{1};

