volScalarField()
fvScalarMatrix SEqn
(
	 fvm::ddt(s)
	+fvm::div(phi,s) 
	-fvm::laplacian(alphaVisc,T)
  
	//-fvm::laplacian(alphaTurb1,T)	
    -fvm::laplacian(alphaTurb2,T)
    
    //-fvm::laplacian(alphaTurb3,T)
    //-fvm::laplacian(alphaTurb4,T)	
	== fvOptions(T)

);

SEqn.relax();
SEqn.solve();
