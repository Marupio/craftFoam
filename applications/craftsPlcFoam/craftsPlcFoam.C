/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    craftsPlcFoam

Description
    Anaerobic Digestion Model with Multi-Dimensional Architecture

    A general coupled reacting flow solver for incompressible flows with user-
    defined implicit routines.  Solves partial differential algebraic equations
    (PDAE) using a point-implicit coupled PDE solver, and user-defined implicit
    routines, which are automatically coupled back into the PDE solver through
    its source term.

SourceFiles
    craftsPlcFoam.C
    createFields.H
    solverFunctions.H
    solverLauncher.H

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "multiSolver.H"
#include "plcEmulator.H"
#include "craftsModel.H"
#include "customUserDefines.H"
#include "solverFunctions.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

    multiSolver multiRun
    (
        multiSolver::multiControlDictName,
        args.rootPath(),
        args.caseName()
    );

    word solverDomain;

    // Find initial solver domain
    timeCluster tcSource(multiRun.initialDataSource());
    
    solverDomain = tcSource.solverDomainName();
    if (solverDomain == "default")
    {
        solverDomain = multiRun.startDomain();
        if (solverDomain == "default")
        {
            FatalIOErrorIn("craftsPlcFoam::main", multiRun.multiControlDict())
                << "Cannot determine start solverDomain. Use keyword "
                << "'startDomain' or change 'initialStartFrom'."
                << exit(FatalIOError);
        }
    }

    // Initialize solver domain
    multiRun.setSolverDomain(solverDomain);

    // Apparently I can't use setSolverDomain twice without running the solver
    bool bugWorkAround(true);

    // Create controller, reads current solverDomain and triggers
    plcEmulator control(multiRun);

    // Main loop (superLoop)
    while (multiRun.run())
    {
        solverDomain = control.currentSolverDomainName();

        if (!bugWorkAround)
        {
            multiRun.setSolverDomain(solverDomain);
        }
        bugWorkAround = false;
        
        if (multiRun.run())
        {
            Info << "*** CRAFTS *** Switching to solverDomain "
                << solverDomain << " ***\n" << endl;

#           undef createPhi_H
#           undef createPhiV_H
#           undef initContinuityErrs_H

#           include "solverLauncher.H"

            multiRun++;
        }
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
