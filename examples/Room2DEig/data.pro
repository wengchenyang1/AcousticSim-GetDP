DefineConstant[
c0 = 343.0,
freqMin = 10,
freqMax = 200,
// Geo:
Ind_Propagation_Domain = 100001,
Ind_Walls = 200003,
// Mesh
N_pt_per_lambda = 15,
Max_lc = 0.3,

NbEigenvalues = 5,
EigenvalShiftRe = (freqMin*2*Pi/c0)^2,
EigenvalShiftIm = 0
];

lambda = c0 / freqMax;
lc = lambda / N_pt_per_lambda;
If(lc > Max_lc)
    lc = Max_lc;
EndIf




