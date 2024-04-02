// GetDP - Copyright (C) 1997-2022 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.

#ifndef F_H
#define F_H

#include <vector>
#include "ProData.h"

/* ------------------------------------------------------------------------ */
/*  Warning: the pointers A and V can be identical. You must                */
/*           use temporary variables in your computations: you can only     */
/*           affect to V at the very last time (when you're sure you will   */
/*           not use A anymore).                                            */
/* ------------------------------------------------------------------------ */

#define F_ARG struct Function *Fct, struct Value *A, struct Value *V

/* F_Analytic */

// using +iwt convention
void F_JFIE_ZPolCyl(F_ARG);
void F_RCS_ZPolCyl(F_ARG);
void F_JFIE_TransZPolCyl(F_ARG);
void F_JFIE_SphTheta(F_ARG);
void F_RCS_SphTheta(F_ARG);
void F_JFIE_SphPhi(F_ARG);
void F_RCS_SphPhi(F_ARG);
void F_CurrentPerfectlyConductingSphere(F_ARG);

// using -iwt convention
void F_ElectricFieldPerfectlyConductingSphereMwt(F_ARG);
void F_ElectricFieldDielectricSphereMwt(F_ARG);
void F_ExactOsrcSolutionPerfectlyConductingSphereMwt(F_ARG);
void F_CurrentPerfectlyConductingSphereMwt(F_ARG);

void F_AcousticFieldSoftSphere(F_ARG);
void F_AcousticFieldSoftSphereABC(F_ARG);
void F_DrAcousticFieldSoftSphere(F_ARG);
void F_RCSSoftSphere(F_ARG);
void F_AcousticFieldHardSphere(F_ARG);
void F_RCSHardSphere(F_ARG);
void F_AcousticFieldSoftCylinder(F_ARG);
void F_AcousticFieldSoftCylinderABC(F_ARG);
void F_DrAcousticFieldSoftCylinder(F_ARG);
void F_RCSSoftCylinder(F_ARG);
void F_AcousticFieldHardCylinder(F_ARG);
void F_AcousticFieldHardCylinderABC(F_ARG);
void F_DthetaAcousticFieldHardCylinder(F_ARG);
void F_RCSHardCylinder(F_ARG);

void F_OSRC_C0(F_ARG);
void F_OSRC_R0(F_ARG);
void F_OSRC_Aj(F_ARG);
void F_OSRC_Bj(F_ARG);

void F_pnm(F_ARG);
void F_unm(F_ARG);
void F_snm(F_ARG);
void F_Xnm(F_ARG);
void F_Ynm(F_ARG);
void F_Znm(F_ARG);
void F_Mnm(F_ARG);
void F_Nnm(F_ARG);

void F_DyadGreenHom(F_ARG);
void F_CurlDyadGreenHom(F_ARG);

/* F_PeWe */

void F_ElastodynamicsCylinderCavity(F_ARG);
void F_ElastodynamicsCylinderWall(F_ARG);
void F_ElastodynamicsCylinderWallS(F_ARG);
void F_ElastodynamicsCylinderWallOut(F_ARG);
void F_ElastodynamicsCylinderWallsOut(F_ARG);
void F_ElastoCylinderWallOutAbc(F_ARG);
void F_ElastoCylinderWallsOutAbc(F_ARG);
void F_ElastoCylinderWallOutAbc2(F_ARG);
void F_ElastoCylinderWallOutAbc2Pade(F_ARG);
void F_ElastoCylinderWallsOutAbc2Pade(F_ARG);

/* F_Geometry */

void F_ProjectPointOnEllipse(F_ARG);
void F_Normal(F_ARG);
void F_NormalSource(F_ARG);
void F_Tangent(F_ARG);
void F_TangentSource(F_ARG);
void F_ElementVol(F_ARG);
void F_SurfaceArea(F_ARG);
void F_GetVolume(F_ARG);
void F_GetNumElement(F_ARG);
void F_GetNumElements(F_ARG);
void F_GetNumNodes(F_ARG);
void F_CellSize(F_ARG);
void F_SquNormEdgeValues(F_ARG);

/* F_Raytracing */

void F_CylinderPhase(F_ARG);
void F_DiamondPhase(F_ARG);

/* F_Math */

void F_Exp(F_ARG);
void F_Log(F_ARG);
void F_Log10(F_ARG);
void F_Sqrt(F_ARG);
void F_Sin(F_ARG);
void F_Asin(F_ARG);
void F_Cos(F_ARG);
void F_Acos(F_ARG);
void F_Tan(F_ARG);
void F_Atan(F_ARG);
void F_Sinh(F_ARG);
void F_Cosh(F_ARG);
void F_Tanh(F_ARG);
void F_Atanh(F_ARG);
void F_Erf(F_ARG);
void F_Fabs(F_ARG);
void F_Abs(F_ARG);
void F_Floor(F_ARG);
void F_Ceil(F_ARG);
void F_Fmod(F_ARG);
void F_Sign(F_ARG);
void F_Min(F_ARG);
void F_Max(F_ARG);
void F_Jn(F_ARG);
void F_JnComplex(F_ARG);
void F_KnComplex(F_ARG);
void F_Yn(F_ARG);
void F_dJn(F_ARG);
void F_dYn(F_ARG);
void F_JnSph(F_ARG);
void F_YnSph(F_ARG);
void F_dJnSph(F_ARG);
void F_dYnSph(F_ARG);

/* F_ExtMath */

void F_Hypot(F_ARG);
void F_Atan2(F_ARG);
void F_TanhC2(F_ARG);
void F_Transpose(F_ARG);
void F_Inv(F_ARG);
void F_Det(F_ARG);
void F_Trace(F_ARG);
void F_RotateXYZ(F_ARG);
void F_Norm(F_ARG);
void F_SquNorm(F_ARG);
void F_Unit(F_ARG);
void F_ScalarUnit(F_ARG);
void F_Cos_wt_p(F_ARG);
void F_Sin_wt_p(F_ARG);
void F_Period(F_ARG);
void F_Interval(F_ARG);
void F_Complex(F_ARG);
void F_Complex_MH(F_ARG);
void F_Re(F_ARG);
void F_Im(F_ARG);
void F_Conj(F_ARG);
void F_Cart2Pol(F_ARG);
void F_Vector(F_ARG);
void F_Tensor(F_ARG);
void F_TensorV(F_ARG);
void F_TensorSym(F_ARG);
void F_TensorDiag(F_ARG);
void F_SquDyadicProduct(F_ARG);
void F_Comp(F_ARG);
void F_CompX(F_ARG);
void F_CompY(F_ARG);
void F_CompZ(F_ARG);
void F_CompXX(F_ARG);
void F_CompXY(F_ARG);
void F_CompXZ(F_ARG);
void F_CompYX(F_ARG);
void F_CompYY(F_ARG);
void F_CompYZ(F_ARG);
void F_CompZX(F_ARG);
void F_CompZY(F_ARG);
void F_CompZZ(F_ARG);
void F_Cart2Sph(F_ARG);
void F_Cart2Cyl(F_ARG);
void F_UnitVectorX(F_ARG);
void F_UnitVectorY(F_ARG);
void F_UnitVectorZ(F_ARG);

/* F_Coord */

// se basent sur le uvw courant (-> en cal)
void F_CoordX(F_ARG);
void F_CoordY(F_ARG);
void F_CoordZ(F_ARG);
void F_CoordXYZ(F_ARG);

// se basent sur le xyz courant, i.e. les coord d'un noeud (-> en pre)
void F_aX_bY_cZ(F_ARG);
void F_aX21_bY21_cZ21(F_ARG);

void F_CoordXS(F_ARG);
void F_CoordYS(F_ARG);
void F_CoordZS(F_ARG);
void F_CoordXYZS(F_ARG);

/* F_Misc */

void F_Printf(F_ARG);
void F_Rand(F_ARG);
void F_CompElementNum(F_ARG);
void F_ElementNum(F_ARG);
void F_QuadraturePointIndex(F_ARG);
void F_GetCpuTime(F_ARG);
void F_GetWallClockTime(F_ARG);
void F_GetMemory(F_ARG);
void F_SetNumberRunTime(F_ARG);
void F_SetNumberRunTimeWithChoices(F_ARG);
void F_GetNumberRunTime(F_ARG);
void F_SetVariable(F_ARG);
void F_SetCumulativeVariable(F_ARG);
void F_GetVariable(F_ARG);
void F_ValueFromTable(F_ARG);
void F_ValueFromFile(F_ARG);
void F_VirtualWork(F_ARG);

void F_Felec(F_ARG);

void F_dFxdux(F_ARG);
void F_dFydux(F_ARG);
void F_dFzdux(F_ARG);
void F_dFxduy(F_ARG);
void F_dFyduy(F_ARG);
void F_dFzduy(F_ARG);
void F_dFxduz(F_ARG);
void F_dFyduz(F_ARG);
void F_dFzduz(F_ARG);

void F_dFxdv(F_ARG);
void F_dFydv(F_ARG);
void F_dFzdv(F_ARG);

void F_dWedxdv(F_ARG);
void F_dWedydv(F_ARG);
void F_dWedzdv(F_ARG);

void F_NodeForceDensity(F_ARG);
void F_AssDiag(F_ARG); /* pour Patrick */

void F_AtIndex(F_ARG);

/* F_Interpolation */

void F_InterpolationLinear(F_ARG);
void F_dInterpolationLinear(F_ARG);
void F_dInterpolationLinear2(F_ARG);
void F_InterpolationAkima(F_ARG);
void F_dInterpolationAkima(F_ARG);
void F_InterpolationBilinear(F_ARG);
void F_dInterpolationBilinear(F_ARG);
void F_InterpolationTrilinear(F_ARG);
void F_dInterpolationTrilinear(F_ARG);
bool Fi_InterpolationBilinear(double *x, double *y, double *M, int NL, int NC,
                              double xp, double yp, double *zp);
bool Fi_dInterpolationBilinear(double *x, double *y, double *M, int NL, int NC,
                               double xp, double yp, double *dzp_dx,
                               double *dzp_dy);
bool Fi_InterpolationTrilinear(double *x, double *y, double *z, double *M,
                               int NX, int NY, int NZ, double xp, double yp,
                               double zp, double *vp);
bool Fi_dInterpolationTrilinear(double *x, double *y, double *z, double *M,
                                int NX, int NY, int NZ, double xp, double yp,
                                double zp, double *dvp_dx, double *dvp_dy,
                                double *dvp_dz);
void Fi_InitListX(F_ARG); // List
void Fi_InitListXY(F_ARG); // ListAlt
void Fi_InitListXY2(F_ARG);
void Fi_InitAkima(F_ARG);
void Fi_InitListMatrix(F_ARG);
void Fi_InitListMatrix3D(F_ARG);
void F_ValueFromIndex(F_ARG);
void F_VectorFromIndex(F_ARG);
void Fi_InitValueFromIndex(F_ARG);
void Fi_InitVectorFromIndex(F_ARG);

void F_TransformTensor(F_ARG); /* pour Tuan */
void F_TransformPerm(F_ARG); /* pour Tuan */
void F_TransformPiezo(F_ARG); /* pour Tuan */
void F_TransformPiezoT(F_ARG); /* pour Tuan */

/* F_Hysteresis */

void F_dhdb_Jiles(F_ARG);
void F_dbdh_Jiles(F_ARG);
void F_h_Jiles(F_ARG);
void F_b_Jiles(F_ARG);

void F_dhdb_Ducharne(F_ARG);
void F_h_Ducharne(F_ARG);
void F_nu_Ducharne(F_ARG);
double Fi_h_Ducharne(double *hi, double *bi, double *M, int NL, int NC,
                     double h0, double b0, double b);

// Energy-Based (EB) Hysteresis Model
// Usefull Mathematical functions:
double Mul_VecVec(const double *v1, const double *v2);
void Mul_TensorVec(const double *M, const double *v, double *Mv,
                   const int transpose_M);
void Mul_TensorSymTensorSym(double *A, double *B, double *C);
void Mul_TensorNonSymTensorNonSym(double *A, double *B, double *C);
void Mul_TensorNonSymTensorSym(double *A, double *B, double *C);
void Mul_TensorSymTensorNonSym(double *A, double *B, double *C);
void Inv_Tensor3x3(double *T, double *invT);
void Inv_TensorSym3x3(double *T, double *invT);

// Anhysteretic curve Characteristics:
double Lang(double nhr, double Ja, double ha);
double dLang(double nhr, double Ja, double ha);
double LangOverx(double nhr, double Ja, double ha);
double dLangOverx(double nhr, double Ja, double ha);
double ILang(double nhr, double Ja, double ha);

double Janhy(double nhr, double Ja, double ha);
double dJanhy(double nhr, double Js, double alpha);
double Xanhy(double nhr, double Js, double alpha);
double dXanhy(double nhr, double Js, double alpha);
double IJanhy(double nhr, double Js, double alpha);
double InvJanhy(double nJ, double Js, double alpha);
double dInvJanhy(double nJ, double Js, double alpha);

double Janhy(double nJ, double Ja, double ha, double Jb, double hb);
double dJanhy(double nJ, double Ja, double ha, double Jb, double hb);
double Xanhy(double nJ, double Ja, double ha, double Jb, double hb);
double dXanhy(double nJ, double Ja, double ha, double Jb, double hb);
double IJanhy(double nJ, double Ja, double ha, double Jb, double hb);
double InvJanhy(double nJ, double Ja, double ha, double Jb, double hb);
double dInvJanhy_hr(double nhr, double Ja, double ha, double Jb, double hb);

double u_hr(double nhr, double Ja, double ha, double Jb, double hb);
double u_J(double nJ, double Js, double alpha);

void Vector_Jk_From_hrk(const double hrk[3], void *params, double Jk[3]);

void Vector_hrk_From_Jk(const double Jk[3], void *params, double hrk[3]);
void Tensor_dJkdhrk(const double hr[3], void *params, double mutg[6]);

// Pseudo-Potential Functional Characteristics:
double fct_omega_VAR(const double h[3], const double Jk[3], void *params);

void fct_d_omega_VAR(const double Jk[3], double *d_omega, double h[3],
                     void *params);

void fct_dd_omega_VAR(const double h[3], const double Jk[3], void *params,
                      double *ddfdJ2);

// Usefull Functions for the Full Differential Approach:
void fct_ehirr_DIFF_3d(const double x[2], double ehi[3]);
void fct_hr_DIFF_3d(const double x[2], const double kappa, const double ehi[3],
                    const double h[3], double xup[3]);
void fct_fall_DIFF_3d(const double ang[2], void *params, double fall[3]);

void fct_hirr_DIFF_2d(double x, double kappa, double ex[3]);
void fct_hr_DIFF_2d(double x, double kappa, double h[3], double xup[3]);
double fct_f_DIFF_2d(double y, void *params);

// Energy-Based Model - Vector Update:
void Vector_Update_Jk_VAR(const double h[3], double Jk[3], double Jk0[3],
                          void *params);

void Vector_Update_hrk_DIFF(const double h[3], double hrk[3], double hrk0[3],
                            void *params);

void Vector_Update_hrk_VPM(const double h[3], double hrk[3], double hrk0[3],
                           void *params);

void Vector_Update_hrk_VPM_(const double h[3], double hr[3],
                            const double hrp[3], const double kappa);

void Vector_b_EB(const double h[3], double b[3], double *Xk_all, void *params);

void Vector_h_EB(const double b[3], double bc[3], double h[3], double *Xk_all,
                 void *params);

// Energy-Based Model - Tensor Construction:
void Tensor_dJkdh_VAR(const double h[3], const double Jk[3], void *params,
                      double *dJkdh);

void Tensor_dJkdh_VPMorDIFF(const double h[3], const double hrk[3],
                            void *params, double *dJkdh);

void Tensor_dhrkdh_DIFF_ana(const double h[3], const double hrk[3],
                            const double dJkdhrk[6], void *params,
                            double *dhrkdh);

void Tensor_dhrkdh_VPM_ana(const double h[3], const double hrk[3], void *params,
                           double *dhrkdh);

void Tensor_dbdh_ana(const double h[3], const double *Xk_all, void *params,
                     double *dbdh);

void Tensor_num(void (*f)(const double *, double *, double *, void *),
                const double x[3], const double delta0, double *Xk_all,
                void *params, double *dfdx);

void Tensor_dhdb_GoodBFGS(const double dx[3], const double df[3], double *dhdb);

// Energy-Based Model - GetDP Functions:
void F_Cell_EB(F_ARG);
void F_h_EB(F_ARG);
void F_b_EB(F_ARG);
void F_hrev_EB(F_ARG);
void F_Jrev_EB(F_ARG);
void F_dbdh_EB(F_ARG);
void F_dhdb_EB(F_ARG);

// Energy-Based Model - GetDP Functions (DEPRECIATED):
void F_Update_Cell_K(F_ARG);
void F_h_Vinch_K(F_ARG);
void F_b_Vinch_K(F_ARG);
void F_hr_Vinch_K(F_ARG);
void F_Jr_Vinch_K(F_ARG);
void F_dbdh_Vinch_K(F_ARG);
void F_dhdb_Vinch_K(F_ARG);

/* F_MultiHar */

void F_MHToTime(F_ARG);

// the following should go somewhere else
void Fi_MHTimeIntegration(int TypePsi, int NbrTimePoint,
                          List_T *WholeQuantity_L, int FreqOffSet,
                          struct Element *Element,
                          struct QuantityStorage *QuantityStorage_P0, double u,
                          double v, double w, struct Value *ValueOut);
void F_MHToTime0(int init, struct Value *A, struct Value *V, int iTime,
                 int NbrTimePoint, double *TimeMH); /* OJO!!! int *init */
void MHTransform(struct Element *Element,
                 struct QuantityStorage *QuantityStorage_P0, double u, double v,
                 double w, std::vector<struct Value> &MH_Inputs,
                 struct Expression *Expression_P, int NbrPoints,
                 struct Value &MH_Output);

/* F_BiotSavart */
void F_BiotSavart(F_ARG);
void F_Pocklington(F_ARG);

/* F_Gmsh */
void F_Field(F_ARG);
void F_ScalarField(F_ARG);
void F_VectorField(F_ARG);
void F_TensorField(F_ARG);
void F_ComplexScalarField(F_ARG);
void F_ComplexVectorField(F_ARG);
void F_ComplexTensorField(F_ARG);
void F_GradScalarField(F_ARG);
void F_GradVectorField(F_ARG);
void F_GradComplexScalarField(F_ARG);
void F_GradComplexVectorField(F_ARG);
void F_Distance(F_ARG);

/* F_DiffGeom */
void F_Hodge(F_ARG);
void F_InnerProduct(F_ARG);
void F_Sharp(F_ARG);
void F_Flat(F_ARG);
void F_WedgeProduct(F_ARG);
void F_TensorProduct(F_ARG);
void F_InteriorProduct(F_ARG);
void F_PullBack(F_ARG);
void F_PullBackMetric(F_ARG);
void F_PushForward(F_ARG);
void F_InvPullBack(F_ARG);
void F_InvPushForward(F_ARG);

/* F_Octave */
void F_Octave(F_ARG);

/* F_Python */
void F_Python(F_ARG);

#endif
