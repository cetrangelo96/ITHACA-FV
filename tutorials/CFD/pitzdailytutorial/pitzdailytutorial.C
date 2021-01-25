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
/*
#ifndef fvCFD_H
#define fvCFD_H

#include "parRun.H" // Helper class for initializing parallel jobs from the command arguments. Also handles cleanup of parallel (or serial) jobs.

#include "Time.H" // Class to control time during OpenFOAM simulations that is also the top-level objectRegistry.

#include "fvMesh.H" // Mesh data needed to do the Finite Volume discretisation.

#include "fvc.H" // Namespace of functions to calculate explicit derivatives.
#include "fvMatrices.H" // A special matrix type and solver, designed for finite volume solutions of scalar equations. 
#include "fvm.H" // Namespace of functions to calculate implicit derivatives returning a matrix.
#include "linear.H" // Central-differencing interpolation scheme class

#include "uniformDimensionedFields.H" // Dimensioned<Type> registered with the database as a registered IOobject which has the functionality of a uniform field and allows values from the top-level code to be passed to boundary conditions etc. Is a 'global' field, same on all processors

#include "calculatedFvPatchFields.H" // This boundary condition is not designed to be evaluated; it is assumed that the value is assigned via field assignment, and not via a call to e.g. \c updateCoeffs or \c evaluate.
#include "extrapolatedCalculatedFvPatchFields.H" // This boundary condition applies a zero-gradient condition from the patch internal field onto the patch faces when \c evaluated but may also be assigned. \c snGrad returns the patch gradient evaluated from the current internal and patch field values rather than returning zero.
#include "fixedValueFvPatchFields.H" // This boundary condition supplies a fixed value constraint, and is the base class for a number of other boundary conditions.
#include "zeroGradientFvPatchFields.H" // This boundary condition applies a zero-gradient condition from the patch internal field onto the patch faces.
#include "fixedFluxPressureFvPatchScalarField.H" //This boundary condition sets the pressure gradient to the provided value such that the flux on the boundary is that specified by the velocity boundary condition.
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H" // For cases which do no have a pressure boundary adjust the balance of fluxes to obey continuity.  Return true if the domain is closed.

#include "findRefCell.H" // Find the reference cell nearest (in index) to the given cell but which is not on a cyclic, symmetry or processor patch.
#include "IOMRFZoneList.H" // List of MRF zones with IO functionality.  MRF zones are specified by a list of dictionary entries

#include "constants.H" // Different types of constants
#include "gravityMeshObject.H" // Gravitational acceleration vector. Although termed a \em MeshObject it is registered on Time only and thus identical for all regions.

#include "columnFvMesh.H" // Generates a 1D column representation of a mesh based on an existing mesh and/or fields

#include "OSspecific.H" // Functions used by OpenFOAM that are specific to POSIX compliant operating systems and need to be replaced or emulated on other systems.
#include "argList.H" // Extract command arguments and options from the supplied \a argc and \a argv parameters.
#include "timeSelector.H" // A List of scalarRange for selecting times. The timeSelector provides a convenient means of selecting multiple times.

#ifndef namespaceFoam
#define namespaceFoam
    using namespace Foam;
#endif

#endif
*/
#include "pisoControl.H" // Specialization of the pimpleControl class for PISO control.

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H" // Execute application functionObjects to post-process existing results.

    #include "addCheckCaseOptions.H" // Check case set-up only using a single time step. Or, Check case set-up and write only using a single time step. 
    #include "setRootCaseLists.H"
    /*
    // This is setRootCase, but with additional solver-related listing
    #include "setRootCaseListOptions.H" // Declare some "standard" list options
    #include "setRootCase.H" // Construct from (int argc, char* argv[]), - use argList::argsMandatory() to decide on checking command arguments. - check validity of the options
    #include "setRootCaseListOutput.H" // Process some "standard" list options
    */
    #include "createTime.H" //Foam::Info<< "Create time\n" << Foam::endl;
		            //Foam::Time runTime(Foam::Time::controlDictName, args);
    #include "createMesh.H" // Create mesh for time = 0

    pisoControl piso(mesh);

    #include "createFields.H"
    /* Reading transportProperties
       Reading field p
       Reading field U
       Reading/calculating face flux field phi*/
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

    return 0;
}

