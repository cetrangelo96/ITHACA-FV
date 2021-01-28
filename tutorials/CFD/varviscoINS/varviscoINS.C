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
    Transient solver for incompressible, laminar flow of Newtonian fluids with 
    variable viscosity.
SourceFiles
    varviscoINS.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H" // Specialization of the pimpleControl class for PISO control.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class varviscoINS
{
    public:
	// Constructor
	explicit varviscoINS(int argc, char* argv[])
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
	    _piso = autoPtr<pisoControl>
              	    (
               	        new pisoControl
                        (
                            mesh
                    	)
              	    );
	#include "createFields.H"
	};
	// Destructor
	~varviscoINS() {};
	/// Reference pressure cell
        label pRefCell;
        /// Reference pressure value
        scalar pRefValue;
	// argList
	autoPtr<argList> _args;
	/// simpleControl
        autoPtr<pisoControl> _piso;
	// Dummy variables to transform icoFoam into a class
	/// Pressure field
        autoPtr<volScalarField> _p;
        /// Velocity field
        autoPtr<volVectorField> _U;
        /// Viscosity
        autoPtr<surfaceScalarField> _nu;
	/// Flux
        autoPtr<surfaceScalarField> _phi;
        /// Mesh
        mutable autoPtr<fvMesh> _mesh;
        /// Time
        autoPtr<Time> _runTime;
	// Functions
	/// Define the viscosity function
        void compute_nu()
        {   fvMesh& mesh = _mesh();
            surfaceScalarField yPos = mesh.Cf().component(vector::Y);
            surfaceScalarField xPos = mesh.Cf().component(vector::X);
	    surfaceScalarField& nu = _nu();
            forAll(xPos, counter)
            {
                 nu[counter] = xPos[counter] + yPos[counter];
            }
        }
        //--------------------------------------------------------------------------
        /// Perform a truthsolve
	void icosolve()
	{
	Time& runTime = _runTime();
        fvMesh& mesh = _mesh();
        volScalarField& p = _p();
        volVectorField& U = _U();
	surfaceScalarField& nu = _nu();
        surfaceScalarField& phi = _phi();
        pisoControl& piso = _piso();
	#include "initContinuityErrs.H"        
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
	
	};
};



int main(int argc, char *argv[])
{
    varviscoINS example( argc, argv);
    example.compute_nu();
    //example.icosolve();
    return 0;
}

