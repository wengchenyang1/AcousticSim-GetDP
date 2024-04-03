Include "data.pro";

// Define the room corners
Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0, 0, lc};
Point(3) = {10, 5, 0, lc};
Point(4) = {10, 8, 0, lc};
Point(5) = {7, 8, 0, lc};
Point(6) = {3, 8, 0, lc};
Point(7) = {0, 8, 0, lc};
Point(8) = {0, 5, 0, lc};

// Connect the corners with lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

// Create a surface from the lines
Line Loop(9) = {1, 2, 3, 4, 5, 6, 7, 8};
ind_surf = 1000;
Plane Surface(ind_surf) = {9};

// Define the propagation domain
Physical Surface(Ind_Propagation_Domain) = {ind_surf};

Physical Curve(Ind_Walls) = {1, 2, 3, 4, 5, 6, 7, 8};

Physical Point(Ind_PrintPoint) = {1};   // Printpoint