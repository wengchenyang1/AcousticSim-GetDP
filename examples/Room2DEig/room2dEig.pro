//========================================================
// Modified based on the following:
//
// GetDP code for simulation of EIGENVALUES PROBLEMS
//   - Helmholtz equation (scalar/vector and 1D/2D/3D)
//   - Nodal/Edge finite-elements (Form0/Form1 for u)
// Contributors: A. Itagi (original version, 2007),
//   B. Kubicek (minor modif), A. Modave (major modif)
//========================================================

Include "data.pro" ;

Group {
  PropDomain = Region[Ind_Propagation_Domain];
  Wall = Region[Ind_Walls];
  PrintPoint = Region[Ind_PrintPoint];
}

Jacobian {
  { Name Jac ;
    Case {
      { Region Wall ; Jacobian Sur ; }
      { Region PropDomain ; Jacobian Vol ; }
      { Region PrintPoint ; Jacobian Sur ; }
    }
  }
}

Integration {
  { Name I1 ;
      Case {
      { Type Gauss ;
          Case {
          { GeoElement Point ; NumberOfPoints  1 ; }
          { GeoElement Line ; NumberOfPoints  4 ; }
          { GeoElement Triangle ; NumberOfPoints  6 ; }
          { GeoElement Quadrangle ; NumberOfPoints 7 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 15 ; }
          { GeoElement Hexahedron ; NumberOfPoints 34 ; }
          }
      }
      }
  }
}

FunctionSpace {
  {Name H_grad; Type Form0;
      BasisFunction{
          {Name u; NameOfCoef ui; Function BF_Node;
          Support Region[{PropDomain, Wall}]; Entity NodesOf[All];}
      }
  }
}

Formulation {
  { Name Helmholtz ; Type FemEquation ;
    Quantity {
      { Name u ; Type Local ; NameOfSpace H_grad ; }
    }
    Equation {
      Galerkin { [ Dof{d u} , {d u} ] ;
                 In PropDomain ; Integration I1 ; Jacobian Jac ; }
      Galerkin { DtDtDof[ Dof{u} , {u} ] ;
                 In PropDomain ; Integration I1 ; Jacobian Jac ; }
    }
  }
}

Function{
  eigfilter = 1e-3;
  EigFilter[] = (Norm[$EigenvalueReal] > eigfilter);
}
Resolution {
   { Name Reso ;
    System {
      { Name A ; NameOfFormulation Helmholtz ;
        Type ComplexValue ;
      }
    }
    Operation {
      CreateDir["output/"] ;
      GenerateSeparate[A] ;
      EigenSolve[A, NbEigenvalues, EigenvalShiftRe, EigenvalShiftIm, EigFilter[]] ;
      SaveSolutions[A] ;
    }
  }
}

PostProcessing {
  { Name PostPro ; NameOfFormulation Helmholtz ;
    Quantity {
      {Name p ; Value{ Local{ [{u}] ; In PropDomain ; Jacobian Jac; } } }
      {Name eigFreq;  Value { Local{ [$EigenvalueReal*c0/2/Pi]; In PrintPoint; Jacobian Jac; } } }
    }
  }
}

PostOperation {
  { Name PostOp ; NameOfPostProcessing PostPro ;
    Operation {
      For n In {0:(NbEigenvalues-1)}
        Print [ p, OnElementsOf PropDomain, TimeStep{n}, File StrCat["output/eigenvector",Sprintf("%g",n),".pos"]];
      EndFor
      Print [eigFreq, OnElementsOf PrintPoint, Format TimeTable, File "output/EigenValuesReal.pos"];
    }
  }
}


DefineConstant[
  P_ = {"PostOp", Name "GetDP/2PostOperationChoices", Visible 0, ReadOnly 1},
  C_ = {"-solve -pos -v2", Name "GetDP/9ComputeCommand", Visible 0}
] ;

