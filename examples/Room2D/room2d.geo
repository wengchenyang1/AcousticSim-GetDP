Include "data.pro";

// Define the room corners
Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0, 0, lc};
Point(3) = {10, 5, 0, lc};
Point(4) = {7, 5, 0, lc};
Point(5) = {7, 8, 0, lc};
Point(6) = {3, 8, 0, lc};
Point(7) = {3, 5, 0, lc};
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

// Define point source
ind_src = 99999;
Point(ind_src) = {X_source, Y_source, 0, lc};
Field[1] = Distance;
Field[1].PointsList = {ind_src};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = Lc_source;
Field[2].SizeMax = lc;
Field[2].DistMin = Lc_source*5;
Field[2].DistMax = Lc_source*90;

// Use the minimum of all the fields as the background mesh size field
Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
