// Lib_Elasticity_u.pro
//
// Template library for elastostatics, elastodynamics and modal analysis using a
// displacement (u) formulation, in both 2D (plane stress or plain strain) and
// 3D.

// Default definitions of constants, groups and functions that can/should be
// redefined from outside the template:

DefineConstant[
  modelPath = "", // default path of the model
  resPath = StrCat[modelPath, "res/"], // path for post-operation files
  modelDim = 2, // default model dimension (2D)
  Flag_PlaneStress = 0, // plain stress in 2D?
  Flag_Regime = 0, // static (0), harmonic (1), time-domain (2), modal (3)
  Freq = 1, // frequency (for harmonic simulations)
  Freq_Target = 1, // frequency target (for modal simulations)
  Num_Modes = 10, // number of modes (for modal simulations)
  TimeInit = 0, // intial time (for time-domain simulations)
  TimeFinal = 1/50, // final time (for time-domain simulations)
  DeltaTime = 1/500, // time step (for time-domain simulations)
  Flag_Axi = 0 // axisymmetric model?
];

Group {
  DefineGroup[
    // Full elastic domain:
    Vol_Mec,

    // Subsets of Vol_Mec:
    Vol_F_Mec, // region with imposed force

    // Boundaries:
    Sur_Neu_Mec // surfaces with Neumann boundary conditions (pressure)
  ];
  Dom_Mec = Region[ {Vol_Mec, Sur_Neu_Mec} ];
}

Function{
  DefineFunction[
    E, // Young modulus (in Vol_Mec)
    nu, // Poisson coefficient (in Vol_Mec)
    rho, // mass density (in Vol_Mec)
    f, // force per unit volume (in Vol_Force_Mec)
    sigman // traction (on Sur_Neu_Mec)
  ];
}

// End of definitions.

Jacobian {
  { Name JacVol_Mec;
    Case {
      If(Flag_Axi && modelDim < 3)
        { Region All; Jacobian VolAxiSqu; }
      Else
        { Region All; Jacobian Vol; }
      EndIf
    }
  }
  { Name JacSur_Mec;
    Case {
      If(Flag_Axi && modelDim < 3)
        { Region All; Jacobian SurAxi; }
      Else
        { Region All; Jacobian Sur; }
      EndIf
    }
  }
}

Integration {
  { Name Int_Mec;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point; NumberOfPoints 1; }
          { GeoElement Line; NumberOfPoints 2; }
          { GeoElement Triangle; NumberOfPoints 3; }
          { GeoElement Quadrangle; NumberOfPoints 4; }
          { GeoElement Tetrahedron; NumberOfPoints 4; }
          { GeoElement Hexahedron; NumberOfPoints 6; }
          { GeoElement Prism; NumberOfPoints 9; }
          { GeoElement Pyramid; NumberOfPoints 8; }

          { GeoElement Line2; NumberOfPoints 3; }
          { GeoElement Triangle2; NumberOfPoints 6; }
          { GeoElement Quadrangle2; NumberOfPoints 7; }
          { GeoElement Tetrahedron2; NumberOfPoints 15; }
          { GeoElement Hexahedron2; NumberOfPoints 34; }
          { GeoElement Prism2; NumberOfPoints 21; }
	}
      }
    }
  }
}

Function {
  If(Flag_PlaneStress) // plane stress (EPC)
    a[] = E[]/(1.-nu[]^2);
    c[] = E[]*nu[]/(1.-nu[]^2);
  Else // plane strain (EPD) or 3D
    a[] = E[]*(1.-nu[])/(1.+nu[])/(1.-2.*nu[]);
    c[] = E[]*nu[]/(1.+nu[])/(1.-2.*nu[]);
  EndIf
  b[] = E[]/2./(1.+nu[]); // = mu = G

  C_xx[] = Tensor[ a[],0  ,0  ,    0  ,b[],0  ,    0  ,0  ,b[] ];
  C_xy[] = Tensor[ 0  ,c[],0  ,    b[],0  ,0  ,    0  ,0  ,0   ];
  C_xz[] = Tensor[ 0  ,0  ,c[],    0  ,0  ,0  ,    b[],0  ,0   ];

  C_yx[] = Tensor[ 0  ,b[],0  ,    c[],0  ,0  ,    0  ,0  ,0   ];
  C_yy[] = Tensor[ b[],0  ,0  ,    0  ,a[],0  ,    0  ,0  ,b[] ];
  C_yz[] = Tensor[ 0  ,0  ,0  ,    0  ,0  ,c[],    0  ,b[],0   ];

  C_zx[] = Tensor[ 0  ,0  ,b[],    0  ,0  ,0  ,    c[],0  ,0   ];
  C_zy[] = Tensor[ 0  ,0  ,0  ,    0  ,0  ,b[],    0  ,c[],0   ];
  C_zz[] = Tensor[ b[],0  ,0  ,    0  ,b[],0  ,    0  ,0  ,a[] ];
}

FunctionSpace {
  { Name H_ux_Mec; Type Form0;
    BasisFunction {
      { Name sxn; NameOfCoef uxn; Function BF_Node;
        Support Dom_Mec; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef uxn; EntityType NodesOf; NameOfConstraint Displacement_x; }
    }
  }
  { Name H_uy_Mec; Type Form0;
    BasisFunction {
      { Name syn; NameOfCoef uyn; Function BF_Node;
        Support Dom_Mec; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef uyn; EntityType NodesOf; NameOfConstraint Displacement_y; }
    }
  }
  { Name H_uz_Mec; Type Form0;
    BasisFunction {
      { Name syn; NameOfCoef uzn; Function BF_Node;
        Support Dom_Mec; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef uzn; EntityType NodesOf; NameOfConstraint Displacement_z; }
    }
  }
}

Formulation {
  { Name Elasticity_u; Type FemEquation;
    Quantity {
      { Name ux; Type Local; NameOfSpace H_ux_Mec; }
      { Name uy; Type Local; NameOfSpace H_uy_Mec; }
      If(modelDim == 3)
        { Name uz; Type Local; NameOfSpace H_uz_Mec; }
      EndIf
    }
    Equation {
      Integral { [ -C_xx[] * Dof{d ux}, {d ux} ];
        In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
      Integral { [ -C_xy[] * Dof{d uy}, {d ux} ];
        In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
      If(modelDim == 3)
        Integral { [ -C_xz[] * Dof{d uz}, {d ux} ];
          In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
      EndIf

      Integral { [ -C_yx[] * Dof{d ux}, {d uy} ];
        In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
      Integral { [ -C_yy[] * Dof{d uy}, {d uy} ];
        In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
      If(modelDim == 3)
        Integral { [ -C_yz[] * Dof{d uz}, {d uy} ];
          In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
      EndIf

      If(modelDim == 3)
        Integral { [ -C_zx[] * Dof{d ux}, {d uz} ];
          In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
        Integral { [ -C_zy[] * Dof{d uy}, {d uz} ];
          In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
        Integral { [ -C_zz[] * Dof{d uz}, {d uz} ];
          In Vol_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
      EndIf

      If(Flag_Regime)
        Integral { DtDtDof [ -rho[] * Dof{ux} , {ux} ];
          In Vol_Mec ; Jacobian JacVol_Mec ; Integration Int_Mec ; }
        Integral { DtDtDof [ -rho[] * Dof{uy} , {uy} ];
          In Vol_Mec ; Jacobian JacVol_Mec ; Integration Int_Mec ; }
        If(modelDim == 3)
          Integral { DtDtDof [ -rho[] * Dof{uz} , {uz} ];
            In Vol_Mec ; Jacobian JacVol_Mec ; Integration Int_Mec ; }
        EndIf
      EndIf

      If(Flag_Regime != 3)
        Integral { [ CompX[f[]] , {ux} ];
          In Vol_F_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
        Integral { [ CompY[f[]] , {uy} ];
          In Vol_F_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
        If(modelDim == 3)
          Integral { [ CompZ[f[]] , {uy} ];
            In Vol_F_Mec; Jacobian JacVol_Mec; Integration Int_Mec; }
        EndIf

        Integral { [ CompX[sigman[]] , {ux} ];
          In Sur_Neu_Mec; Jacobian JacSur_Mec; Integration Int_Mec; }
        Integral { [ CompY[sigman[]] , {uy} ];
          In Sur_Neu_Mec; Jacobian JacSur_Mec; Integration Int_Mec; }
        If(modelDim == 3)
          Integral { [ CompZ[sigman[]] , {uz} ];
            In Sur_Neu_Mec; Jacobian JacSur_Mec; Integration Int_Mec; }
        EndIf
      EndIf
    }
  }
}

Resolution {
  { Name Elasticity_u;
    System {
      { Name A; NameOfFormulation Elasticity_u;
        If(Flag_Regime == 1)
          Type Complex; Frequency Freq;
        EndIf
      }
    }
    Operation {
      // faster assembly with 2nd order elements
      SetGlobalSolverOptions["-petsc_prealloc 500"];

      If(Flag_Regime == 0 || Flag_Regime == 1)
        Generate[A]; Solve[A]; SaveSolution[A];
      ElseIf(Flag_Regime == 2)
        InitSolution[A]; InitSolution[A] ;
        TimeLoopNewmark[TimeInit, TimeFinal, DeltaTime, 1/4, 1/2] {
          Generate[A]; Solve[A]; SaveSolution[A];
        }
      Else
        GenerateSeparate[A]; EigenSolve[A, Num_Modes, (2*Pi*Freq_Target)^2, 0];
        SaveSolutions[A];
      EndIf
    }
  }
}

PostProcessing {
  { Name Elasticity_u; NameOfFormulation Elasticity_u;
    PostQuantity {
      { Name u; Value {
          If(modelDim == 3)
            Term { [ Vector[ {ux}, {uy}, {uz} ]];
              In Vol_Mec; Jacobian JacVol_Mec; }
          Else
            Term { [ Vector[ {ux}, {uy}, 0 ]];
              In Vol_Mec; Jacobian JacVol_Mec; }
          EndIf
        }
      }
      { Name sigma; Value {
          If(modelDim == 3)
            Term { [ TensorV[ C_xx[]*{d ux} + C_xy[]*{d uy} + C_xz[]*{d uz},
                              C_yx[]*{d ux} + C_yy[]*{d uy} + C_yz[]*{d uz},
                              C_zx[]*{d ux} + C_zy[]*{d uy} + C_zz[]*{d uz} ] ];
              In Vol_Mec; Jacobian JacVol_Mec; }
          Else
            Term { [ TensorV[ C_xx[]*{d ux} + C_xy[]*{d uy},
                              C_yx[]*{d ux} + C_yy[]*{d uy},
                              Vector[0,0,0]] ];
              In Vol_Mec; Jacobian JacVol_Mec; }
          EndIf
          }
      }
    }
  }
}

PostOperation {
  { Name Elasticity_u; NameOfPostProcessing Elasticity_u;
    Operation {
      CreateDir[resPath];
      Print[ sigma, OnElementsOf Vol_Mec, File StrCat[resPath, "sigma.pos"] ];
      Print[ u, OnElementsOf Vol_Mec, File StrCat[resPath, "u.pos"] ];
    }
  }
}
