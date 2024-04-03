/**
See example at
https://onelab.info/onelab_wiki/index.php?title=Tutorial/Laplace_equation_with_Neumann_boundary_condition&mobileaction=toggle_view_desktop
https://www.personal.reading.ac.uk/~sms03snc/fem_notes.pdf
 */

Include "data.pro";

k = 2*Pi*freq/c0;

// =======
// GROUPS
// =======
Group{
    PropDomain = Region[Ind_Propagation_Domain];
    Wall = Region[Ind_Walls];
}

// =========
// JACOBIAN
// =========
Jacobian {
    { Name JVol ; Case { { Region All ; Jacobian Vol ; } } }
    { Name JSur ; Case { { Region All ; Jacobian Sur ; } } }
}
  
  // ======================
  // INTEGRATION PARAMETERS
  // ======================
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

Constraint{
    //Dirichlet boundary condition
    { Name Dirichlet; Type Assign; Case{
        { Region Wall; Value 0.0; }
    } }
}

// Function space
FunctionSpace {
{Name H_grad; Type Form0;
    BasisFunction{
        {Name u; NameOfCoef ui; Function BF_Node;
        Support Region[{PropDomain, Wall}]; Entity NodesOf[All];}
    }
    // Constraint{
    //     //Dirichlet boundary condition (sound-soft)
    //     {NameOfCoef ui; EntityType NodesOf;
    //     NameOfConstraint Dirichlet;}
    // }
}
}

Function {
    SourceAmplitude = 1;

    Sigma = Lc_source*5;
    R_squared[] = (X[]-X_source)*(X[]-X_source) + (Y[]-Y_source)*(Y[]-Y_source);
    SoundSource[] = SourceAmplitude * Exp[-R_squared[] / (2 * Sigma^2)];

    // Robin boundary condition: dp/dn + beta*p = g
    I[] = Complex[0., 1.]; // sqrt(-1)
    Robin_beta[] = I[] * k * 0.0;  // Set to zero to have Neumann boundary conditioin
    Robin_g[] = 0.0;

    Kc[] = Complex[k, -damping];
    K2[] = Kc[] * Kc[];
}

// ============
// FORMULATIONS
// Dof{u}: "Degree Of Freedom". This is used to specify that the quantity is the unknown.
// If "Dof" is not written, then "u" is seen as a test function and not as the unknown.
// ============
Formulation {
{
Name Helmholtz;
Type FemEquation;
Quantity {
    {Name u; Type Local; NameOfSpace H_grad;}
}
Equation {
    //Helmholtz equation
    Galerkin{ [Dof{Grad u}, {Grad u}];
	    In PropDomain; Jacobian JVol; Integration I1;}
    Galerkin{ [-K2[]*Dof{u}, {u}];
	    In PropDomain; Jacobian JVol; Integration I1;}
    Galerkin{ [-SoundSource[], {u}];
        In PropDomain; Jacobian JVol; Integration I1;}
    Galerkin {[Robin_beta[]*Dof{u}, {u}];
        In Wall; Jacobian JSur; Integration I1;}
    Galerkin {[Robin_g[], {u}];
        In Wall; Jacobian JSur; Integration I1;}
}
}
}

// ===========
// RESOLUTIONS
// ===========
Resolution{
{Name Room2D;
    System{ {Name A; NameOfFormulation Helmholtz; Type Complex; } }
    Operation{
    Generate[A]; Solve[A];
    }
}
} // End Resolution



// ================
// POST-PROCESSINGS
// ================
PostProcessing{
{Name Wave; NameOfFormulation Helmholtz;
    Quantity {
    {Name p; Value {Local { [{u}] ; In PropDomain; Jacobian JVol; }}}
    {Name pNorm; Value {Local { [Norm[{u}]] ; In PropDomain; Jacobian JVol; }}}
    }
}
} // End Postprocessing.
  
// ===============
// POST-OPERATIONS
// ===============

PostOperation{
{Name Wave; NameOfPostProcessing Wave ;
Operation {
    Print [p, OnElementsOf PropDomain, File "u.pos"];
    Print [pNorm, OnElementsOf PropDomain, File "u_Abs.pos"];
}
}
} // End PostOperation
