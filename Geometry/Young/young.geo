// Square obstacle
// List of points
d = 0.1;
x = 1;
y = 1;
lSlit=0.01;
eSlit=0.1;
ySlit=0.1;

Point(1) = {0, 0, 0, d};
Point(2) = {x, 0, 0, d};
Point(3) = {0, y, 0, d};
Point(4) = {x, y, 0, d};

Point(5) = {x/2+lSlit, 0, 0, d};
Point(6) = {x/2-lSlit, 0, 0, d};
Point(7) = {x/2+lSlit, y, 0, d};
Point(8) = {x/2-lSlit, y, 0, d};

Point(13) = {x/2+lSlit, 0+ySlit+eSlit, 0, d/2};
Point(14) = {x/2-lSlit, 0+ySlit+eSlit, 0, d/2};
Point(15) = {x/2+lSlit, y-ySlit-eSlit, 0, d/2};
Point(16) = {x/2-lSlit, y-ySlit-eSlit, 0, d/2};

Point(9) = {x/2+lSlit, 0+ySlit, 0, d/2};
Point(10) = {x/2-lSlit, 0+ySlit, 0, d/2};
Point(11) = {x/2+lSlit, y-ySlit, 0, d/2};
Point(12) = {x/2-lSlit, y-ySlit, 0, d/2};


// List of lines
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 8};
//+
Line(3) = {8, 12};
//+
Line(4) = {12, 11};
//+
Line(5) = {11, 7};
//+
Line(6) = {7, 4};
//+
Line(7) = {4, 2};
//+
Line(8) = {2, 5};
//+
Line(9) = {5, 9};
//+
Line(10) = {9, 10};
//+
Line(11) = {10, 6};
//+
Line(12) = {6, 1};
//+
Line(13) = {16, 14};
//+
Line(14) = {14, 13};
//+
Line(15) = {13, 15};
//+
Line(16) = {15, 16};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Curve Loop(2) = {13, 14, 15, 16};

Plane Surface(1) = {1, 2};
Physical Surface(1) = {1};

// Boundary conditions
Physical Curve("BC_Left") = {1};
Physical Curve("BC_Up_1") = {2};
Physical Curve("BC_Up_2") = {6};
Physical Curve("BC_Down_1") = {12};
Physical Curve("BC_Down_2") = {8};
Physical Curve("BC_Right") = {7};

Physical Curve("BC_Obstacle_Up") = {3,4,5};
Physical Curve("BC_Obstacle") = {13, 14, 15 ,16};
Physical Curve("BC_Obstacle_Down") = {9,10,11};

// Transfinite Surface{1};