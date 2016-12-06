Point(1) = {0, 0, 0};
Point(2) = {0, 5, 0};
Point(3) = {1, 5, 0};
Point(4) = {1, 0, 0};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};

Plane Surface(6) = {5};

Physical Line("bottom") = {4};
Physical Line("top") = {2};
Physical Line("left") = {1};
Physical Line("right") = {3};

Physical Surface(7) = {6};