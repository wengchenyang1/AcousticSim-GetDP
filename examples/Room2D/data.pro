MENU_INPUT = "Input";
DefineConstant[
freq = {10, Min 10, Step 5, Max 100,
    Name StrCat[MENU_INPUT, "/1Frequency"]}

c0 = 343.0,
damping = 0.0006, // At 1000Hz
// Source location
Y_source = 4.0,
X_source = 5.0,
// Geo:
Ind_Propagation_Domain = 100001,
Ind_Walls = 200003,
// Mesh
Lc_source = 0.01,
N_pt_per_lambda = 15,
Max_lc = 0.3
];

lambda = c0 / freq;
lc = lambda / N_pt_per_lambda;
If(lc > Max_lc)
    lc = Max_lc;
EndIf
If(lc < Lc_source)
    lc = Lc_source;
EndIf




