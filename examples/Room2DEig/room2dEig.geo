Include "data.pro";

// Define the room corners
Point(1) = {2.469172954559326, -2.6659555435180664, 0, lc};
Point(2) = {-0.36817723512649536, -3.5220680236816406, 0, lc};
Point(3) = {-1.0533497333526611, -1.25125253200531, 0, lc};
Point(4) = {-3.396975517272949, -1.958394169807434, 0, lc};
Point(5) = {-4.440656661987305, 1.5005990266799927, 0, lc};
Point(6) = {0.7403196096420288, 3.0638532638549805, 0, lc};

// Connect the corners with lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// Create a surface from the lines
Line Loop(9) = {1, 2, 3, 4, 5, 6};
ind_surf = 1000;
Plane Surface(ind_surf) = {9};

// Define the propagation domain
Physical Surface(Ind_Propagation_Domain) = {ind_surf};

Physical Curve(Ind_Walls) = {1, 2, 3, 4, 5, 6};

Physical Point(Ind_PrintPoint) = {1};   // Printpoint//+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Show "*";
