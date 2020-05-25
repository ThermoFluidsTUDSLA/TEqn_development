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
forAll(R_diag,i)
{
    R_diag[i].xx() = max(R[i].xx(),1e-10);
    R_diag[i].yy() = max(R[i].yy(),1e-10);
    R_diag[i].zz() = max(R[i].zz(),1e-10);
}

// Optionally Limit off-diagonal Reynolds stresses according to the the realizibility constraints
// ----------------------------------------------------------------------------------------------
volScalarField maxMag_R_xy("maxMag_R_xy",Foam::sqrt(mag(R_diag.component(symmTensor::XX)*R_diag.component(symmTensor::YY))));
volScalarField maxMag_R_xz("maxMag_R_xz",Foam::sqrt(mag(R_diag.component(symmTensor::XX)*R_diag.component(symmTensor::ZZ))));
volScalarField maxMag_R_yz("maxMag_R_yz",Foam::sqrt(mag(R_diag.component(symmTensor::YY)*R_diag.component(symmTensor::ZZ))));

forAll(R_dev,i)
{
        if(std::signbit(R_dev[i].xy())) // Negative value
        {
            R_dev[i].xy() = max(R_dev[i].xy(),-1*maxMag_R_xy[i]);
        }
        else
        {
            R_dev[i].xy() = min(R_dev[i].xy(),+1*maxMag_R_xy[i]);
        }
        if(std::signbit(R[i].xz())) // Negative value
        {
            R_dev[i].xz() = max(R_dev[i].xz(),-1*maxMag_R_xz[i]);
        }
        else
        {
            R_dev[i].xz() = min(R_dev[i].xz(),+1*maxMag_R_xz[i]);
        }
        if(std::signbit(R[i].xz())) // Negative value
        {
            R_dev[i].yz() = max(R_dev[i].yz(),-1*maxMag_R_yz[i]);
        }
        else
        {
            R_dev[i].yz() = min(R_dev[i].yz(),+1*maxMag_R_yz[i]);
        }
}

volSymmTensorField ReynoldsStresses = R_diag + R_dev;

// Turbulent (subgrid) production 
// ------------------------------
//volSymmTensorField P(-twoSymm(R & fvc::grad(U)));

// Calculation of turbulent-diffusivitiy tensors ::
// ------------------------------------------------

dimensionedSymmTensor I("I",dimensionSet(0,0,0,0,0,0,0),symmTensor(1,0,0,1,0,1));

volSymmTensorField alphaTurb1(
                              "alphaTurb1",
                              turbulence->nut()*I/0.9
                             );

volSymmTensorField alphaTurb2(
                              "alphaTurb2",
                              timeScale*(
                                          C2*ReynoldsStresses//+
                                        )
                              );
/*
volTensorField  alphaTurb3(
                              "alphaTurb3",
                              timeScale*(
                                          C3*timeScale*boundedK*fvc::grad(U)
                                        )
                          );

volSymmTensorField  alphaTurb4(
                              "alphaTurb4",
                              timeScale*(
                                          C4*timeScale*P
                                        )
                          );

*/

// Solving the temperature equation ::
// ***********************************

fvScalarMatrix TEqn
(
	fvm::ddt(T)
	+fvm::div(phi,T) 
	-fvm::laplacian(alphaVisc,T)
	//-fvm::laplacian(alphaTurb1,T)	
    -fvm::laplacian(alphaTurb2,T)
    //-fvm::laplacian(alphaTurb3,T)
    //-fvm::laplacian(alphaTurb4,T)	
	== fvOptions(T)

);

TEqn.relax();
TEqn.solve();

Info << "***********************************************" << endl;
Info << "Maximum ::" << gMax(T) << endl;
Info << "Minimum ::" << gMin(T) << endl;
Info << "***********************************************" << endl;









