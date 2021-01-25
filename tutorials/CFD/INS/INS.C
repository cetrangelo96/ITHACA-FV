/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Description
    icoFoam = Transient solver for incompressible, laminar flow of Newtonian fluids.
SourceFiles
    pitzdailytutorial.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H" // Specialization of the pimpleControl class for PISO control.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class INS
{
    public:
	// Constructor
	explicit INS(int argc, char* argv[])
	    /*:
	    U(),
	    p(),
            nu()*/
	{
	    _args = autoPtr<argList>
                (
                    new argList(argc, argv)
                );

            if (!_args->checkRootCase())
            {
                Foam::FatalError.exit();
            }

            argList& args = _args();
        #include "createTime.H"
        #include "createMesh.H"
	    // pisoControl piso(mesh);
	#include "createFields.H"
	/*argList::addNote
       (
           "Transient solver for incompressible, laminar flow"
           " of Newtonian fluids."
       ); 

       #include "postProcess.H" // Execute application functionObjects to post-process existing results.
       #include "addCheckCaseOptions.H" // Check case set-up only using a single time step. Or, Check case set-up and write only using a single time step. 
       #include "setRootCaseLists.H"
    
       #include "createTime.H" 
       #include "createMesh.H" // Create mesh for time = 0

       pisoControl piso(mesh);

       #include "createFields.H"
    
       #include "initContinuityErrs.H" // Declare and initialise the cumulative continuity error. */
	};
	// Destructor
	~INS() {};
	// argList
	autoPtr<argList> _args;
	// Dummy variables to transform icoFoam into a class
        /// Velocity field
        volVectorField U;
        /// Pressure field
        volScalarField p;
        /// Viscosity
        dimensionedScalar nu;
        /// Mesh
        mutable autoPtr<fvMesh> mesh;
        /// Time
        autoPtr<Time> runTime;
	// Functions

        //--------------------------------------------------------------------------
        /// Perform a truthsolve
	/*void icosolve()
	{
        Info<< "\nStarting time loop\n" << endl;

        while (runTime.loop())
            {
        	Info<< "Time = " << runTime.timeName() << nl << endl;

        	#include "CourantNo.H" // Calculates and outputs the mean and maximum Courant Numbers.

        	// Momentum predictor

        	fvVectorMatrix UEqn
        	(
        	    fvm::ddt(U)
        	  + fvm::div(phi, U)
        	  - fvm::laplacian(nu, U)
        	);

	        if (piso.momentumPredictor())
	        {
	            solve(UEqn == -fvc::grad(p));
	        }
	
	        // --- PISO loop
	        while (piso.correct())
	        {
	            volScalarField rAU(1.0/UEqn.A());
	            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
	            surfaceScalarField phiHbyA
	            (
	                "phiHbyA",
	                fvc::flux(HbyA)
	              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
	            );
	
	            adjustPhi(phiHbyA, U, p);
	
	            // Update the pressure BCs to ensure flux consistency
	            constrainPressure(p, U, phiHbyA, rAU);
	
	            // Non-orthogonal pressure corrector loop
	            while (piso.correctNonOrthogonal())
	            {
	                // Pressure corrector
	
	                fvScalarMatrix pEqn
	                (
	                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
	                );
	
	                pEqn.setReference(pRefCell, pRefValue);
	
	                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
	
	                if (piso.finalNonOrthogonalIter())
	                {
	                    phi = phiHbyA - pEqn.flux();
	                }
	            }
	
	            #include "continuityErrs.H"
	
	            U = HbyA - rAU*fvc::grad(p);
	            U.correctBoundaryConditions();
	        }
	
	        runTime.write();
	
	        runTime.printExecutionTime(Info);
	    }
	
	    Info<< "End\n" << endl;
	
	};*/
};



int main(int argc, char *argv[])
{/*
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H" // Execute application functionObjects to post-process existing results.

    #include "addCheckCaseOptions.H" // Check case set-up only using a single time step. Or, Check case set-up and write only using a single time step. 
    #include "setRootCaseLists.H"
    
    #include "createTime.H" 
    #include "createMesh.H" // Create mesh for time = 0

    pisoControl piso(mesh);

    #include "createFields.H"
    
    #include "initContinuityErrs.H" // Declare and initialise the cumulative continuity error.

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H" // Calculates and outputs the mean and maximum Courant Numbers.

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;*/
}

