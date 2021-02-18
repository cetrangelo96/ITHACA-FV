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
    varviscosteadyNS.C
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class varviscosteadyINS
{
    public:
	// Constructor
	explicit varviscosteadyINS(int argc, char* argv[])
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
	    _simple = autoPtr<simpleControl>
	          (
		      new simpleControl
		      (
		          mesh
		      )
	          );
	    simpleControl& simple = _simple();
	#include "createFields.H"
	#include "createFvOptions.H"
	    turbulence->validate();
	};
	// Destructor
	~varviscosteadyINS() {};
	/// Reference pressure cell
        label pRefCell;
        /// Reference pressure value
        scalar pRefValue;
	// argList
	autoPtr<argList> _args;
	/// simpleControl
        autoPtr<simpleControl> _simple;
	/// Turbulence model
        autoPtr<incompressible::turbulenceModel> turbulence;
        /// Laminar transport (used by turbulence model)
        autoPtr<singlePhaseTransportModel> _laminarTransport;
        /// MRF variable
        autoPtr<IOMRFZoneList> _MRF;
	// Dummy variables to transform icoFoam into a class
	/// Pressure field
        autoPtr<volScalarField> _p;
        /// Velocity field
        autoPtr<volVectorField> _U;
	/// Flux
        autoPtr<surfaceScalarField> _phi;
        /// Mesh
        mutable autoPtr<fvMesh> _mesh;
	/// fvOptions
        autoPtr<fv::options> _fvOptions;
        /// Time
        autoPtr<Time> _runTime;
	// Functions
	/// Define the viscosity function
        void calc_nu()
        {   fvMesh& mesh = _mesh();
            volScalarField yPos = mesh.C().component(vector::Y);
            volScalarField xPos = mesh.C().component(vector::X);
	    const volScalarField& nu = _laminarTransport().nu();
	    volScalarField& nu_ = const_cast<volScalarField&>(nu);
	    forAll(xPos, counter)
            {
		 nu_[counter] = (1 + 6*pow(xPos[counter],2)+xPos[counter]/(1+2*pow(yPos[counter],2)));		 
            }
	    for(label j=0; j<nu_.boundaryField().size(); j++)
            {
		for (label i = 0; i < nu_.boundaryField()[j].size(); i++)
		    {
            	    nu_.boundaryFieldRef()[j][i] = (1 + 6*pow(xPos.boundaryField()[j][i],2)+xPos.boundaryField()[j][i]/(1+2*pow(yPos.boundaryField()[j][i],2)));
		    }
	    }	

        }
        //--------------------------------------------------------------------------
        /// Perform a SIMPLE solver
	void SIMPLEsolve()
	{
	    Time& runTime = _runTime();
            fvMesh& mesh = _mesh();
            volScalarField& p = _p();
            volVectorField& U = _U();
            surfaceScalarField& phi = _phi();
            fv::options& fvOptions = _fvOptions();
            simpleControl& simple = _simple();
            IOMRFZoneList& MRF = _MRF();
    	    singlePhaseTransportModel& laminarTransport = _laminarTransport();
	     #include "initContinuityErrs.H"

	    turbulence->validate();

	    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	    Info<< "\nStarting time loop\n" << endl;

	    while (simple.loop())
	    {
		Info<< "Time = " << runTime.timeName() << nl << endl;

		// --- Pressure-velocity SIMPLE corrector
		{
		    #include "UEqn.H"
		    #include "pEqn.H"
		}

		laminarTransport.correct();
		turbulence->correct();

		runTime.write();

		runTime.printExecutionTime(Info);
	    }

	    Info<< "End\n" << endl;

	};
};



int main(int argc, char *argv[])
{
    varviscosteadyINS example( argc, argv);
    example.calc_nu();
    example.SIMPLEsolve();
    return 0;
}

