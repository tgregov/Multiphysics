// Square obstacle
// List of points
d = 0.03;
x = 1;
y = 1;
xA = 0.45;
yA = 0.5;
xB = 0.55;
yB = 0.4;
xC = 0.55;
yC = 0.6;

Point(1) = {0, 0, 0, d};
Point(2) = {x, 0, 0, d};
Point(3) = {0, y, 0, d};
Point(4) = {x, y, 0, d};
Point(5) = {xA, yA, 0, d};
Point(6) = {xB, yB, 0, d};
Point(7) = {xC, yC, 0, d};


// List of lines
Line(1) = {3, 1};
Line(2) = {4, 3};
Line(3) = {2, 4};
Line(4) = {1, 2};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 5};

// Surface
Curve Loop(1) = {1, 4, 3, 2};
Curve Loop(2) = {5, 6, 7};
Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};

// Boundary conditions
Physical Curve("BC_Left") = {1};
Physical Curve("BC_Up") = {2};
Physical Curve("BC_Down") = {4};
Physical Curve("BC_Right") = {3};

Physical Curve("BC_Obstacle") = {5, 6, 7};

// Transfinite Surface{1};

