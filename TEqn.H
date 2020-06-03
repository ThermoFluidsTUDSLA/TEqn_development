// Additionally add following lines to createFields.H
/* 

// Addition for the temperature equation
// *************************************************************************

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field for the entropy production\n" << endl;
volScalarField Sigma_alpha
(
    IOobject
    (
        "Sigma_alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,2,-3,-1,0,0,0)
);


Info<< "Reading the model for subGridFluxes\n" << endl;

IOdictionary turbulencePropertiesDict
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensioned<scalar> alpha_t_Coeff
(
	dimensioned<scalar>::lookupOrAddToDict
	(
	   "alpha_t_Coeff",
	   turbulencePropertiesDict.subDict("subGridFluxes"),
	   1
	)
);

dimensioned<scalar> Prt
(
	dimensioned<scalar>::lookupOrAddToDict
	(
	   "Prt",
	   turbulencePropertiesDict.subDict("subGridFluxes"),
	   0.9
	)
);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportPropertiesDict
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Pr("Pr",dimensionSet(0,0,0,0,0,0,0),transportPropertiesDict); // Molecular Prandtl number

dimensionedScalar cp("cp",dimensionSet(0,2,-2,-1,0,0),transportPropertiesDict); // Specific isobaric heat capacity

// *************************************************************************

*/


// Solve Temperature equation

// ******************************************************************************************************************************
// ******************************************************************************************************************************

// Calculation of kinematic molecular thermal diffusivity [m2/s]

// ******************************************************************************************************************************
// ******************************************************************************************************************************

volScalarField alphaVisc = turbulence->nu()/Pr;

// ******************************************************************************************************************************
// ******************************************************************************************************************************

// Calculation of turbulent (subgrid) kinematic thermal diffusivity [m2/s]

// ******************************************************************************************************************************
// ******************************************************************************************************************************
//dimensionedScalar C1("C1",dimensionSet(0,0,0,0,0,0,0),-4.55e-2);

dimensionedScalar C2("C2",dimensionSet(0,0,0,0,0,0,0),+3.73e-1);

//dimensionedScalar C3("C3",dimensionSet(0,0,0,0,0,0,0),-3.73e-3);

//dimensionedScalar C4("C4",dimensionSet(0,0,0,0,0,0,0),-2.35e-2);


// Load turbulent quantities 
// -------------------------
dimensionedScalar kMin("kMin",dimensionSet(0,2,-2,0,0,0,0),1e-10);

dimensionedScalar epsilonMin("epsilonMin",dimensionSet(0,2,-3,0,0,0,0),1e-10); 

volScalarField boundedK("boundedK",max(turbulence->k(),kMin));

volScalarField boundedEpsilon(
                               "boundedEpsilon",
                                max(
                                    turbulence->epsilon()+0.5*fvc::laplacian(turbulence->nu(),boundedK),
                                    epsilonMin
                                   )
                             );

volScalarField timeScale(
                          "timeScale",
                           max(
                               boundedK/boundedEpsilon,
                               Foam::sqrt(turbulence->nu()/boundedEpsilon) // Kolmogorov time-scale
                              )
                        );


// Reynolds stresses 
// -----------------

volSymmTensorField R(turbulence->R());

volSymmTensorField R_diag = 0*R;
volSymmTensorField R_dev  = 0*R;

forAll(R,i)
{
    R_diag[i].xx() = R[i].xx();
    R_diag[i].yy() = R[i].yy();
    R_diag[i].zz() = R[i].zz();

    R_dev[i].xy() = R[i].xy();
    R_dev[i].xz() = R[i].xz();
    R_dev[i].yz() = R[i].yz();
}

// Limit normal Reynolds stresses according to the the realizibility constraints
// -----------------------------------------------------------------------------
R_diag.max
(
    dimensionedSymmTensor
    (
        "zero",
        R.dimensions(),
        symmTensor
        (
            1e-10, 0.0, 0.0,
            1e-10, 0.0,
            1e-10
        )
    )
);

volSymmTensorField ReynoldsStresses = R_diag + R_dev;

// Turbulent (subgrid) production 
// ------------------------------
//volSymmTensorField P(-twoSymm(R & fvc::grad(U)));

// Calculation of turbulent-diffusivitiy tensors ::
// ------------------------------------------------

dimensionedSymmTensor I("I",dimensionSet(0,0,0,0,0,0,0),symmTensor(1,0,0,1,0,1));

volSymmTensorField alphaTurbSGH(
                              "alphaTurb1",
                              turbulence->nut()*I/Prt.value()
                             );

volSymmTensorField alphaTurbGGH(
                              "alphaTurb2",
                              C2*timeScale*ReynoldsStresses          
                              );

volSymmTensorField alphaTurb("alphaTurb",alphaTurbSGH);

if(alpha_t_Coeff.value() == 1.0)
{
  Info << "------------------------------------------" << endl;
  Info << "Using SGH hypotesys for the subGrid fluxes" << endl;
  Info << "------------------------------------------" << endl;
}
else
{
  Info << "------------------------------------------" << endl;
  Info << "Using GGH hypotesys for the subGrid fluxes" << endl;
  Info << "------------------------------------------" << endl; 

  alphaTurb = alphaTurbGGH;
}

// Solving the temperature equation ::
// ***********************************


fvScalarMatrix TEqn
(
  fvm::ddt(T)
  +fvm::div(phi,T) 
  -fvm::laplacian(alphaVisc,T)
  -fvm::laplacian(alphaTurb,T) 
  == fvOptions(T)

);

TEqn.relax();
TEqn.solve();

Info << "***********************************************" << endl;
Info << "Maximum T::" << gMax(T) << endl;
Info << "Minimum T::" << gMin(T) << endl;
Info << "***********************************************" << endl;

// Additional check for the entropy production due to the temperature diffusion (laminar + turbulent)

Sigma_alpha = ( alphaVisc*cp/sqr(T) )*( (fvc::grad(T) &fvc::grad(T)) + ( (alphaTurb & fvc::grad(T)) & fvc::grad(T) )/alphaVisc);

scalar min_Sigma_alpha = gMin(Sigma_alpha);

Info << "------------------------------------------" << endl;
Info << "Minimal production of entropy :: " << min_Sigma_alpha << endl;
Info << "------------------------------------------" << endl;

if(std::signbit(min_Sigma_alpha))
{
  Info << "Negative production of entropy, writing the field!!" << endl;
  Sigma_alpha.write();
}






