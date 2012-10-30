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

\*---------------------------------------------------------------------------*/

int findTheMatrixSize(Time& runTime)
{
    int returnMe(0);
    IOdictionary variableDict
    (
        IOobject
        (
            "admVariableDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    wordList varNames(variableDict.toc());
    forAll(varNames, nameI)
    {
        if (varNames[nameI] == "defaults") continue;
        word varType
        (
            variableDict.subDict(varNames[nameI]).lookup("type")
        );
        
        if (varType == "standard" || varType == "implicit")
        {
            returnMe++;
        }
    }
    return returnMe;
}


template<int matrixSize>
void solve
(
    admTime& runTime,
    fvMesh& mesh,
    plcEmulator& control
)
{
#   include "createFields.H"

    // Allow control to adjust the timestep if necessary, also, checks for
    // changes in solverDomain - this will intercept when a run is resumed
    // and switch to the correct solverDomain.
    control.initializeAndAdjustTime(runTime);

    // Timestep control:
        // Solve one double-sized step
        // Save internal fields only
        // Solve two single-sized steps
        // Residual compare
            // If residual > setting, reduce timestep, return to saved state
            // If residual < setting, increase timestep, play on

    while (runTime.run())
    {
        // Save state so we can return here for the two half-timesteps.
        model.saveState(0);

        // True when difference between coarse step and fine double step is
        // acceptable
        bool solutionAcceptable(false);
        
        // Becomes false only when we successfully perform all steps (a coarse
        // step and two fine steps)
        bool badStep(true);
        
        // For adaptive timestepping
        scalar nextDeltaT(runTime.deltaT().value());
        scalar oldNextDeltaT(nextDeltaT);

        // Update error scaling factors
        model.calculateScales();

        while (!solutionAcceptable)
        {
            do // 'loop' executes only once - used to jump to end upon failure
            {
                // Double the timestep, and increment time value
                runTime.doubleDeltaT();
                
                Info << "\n---CRAFTS--- Coarse step"
                    << endl;
                // Solve a single step size
                if (model.coarseStep() != 0)
                {
                    break;
                }
                
                // Return to the original starting point
                model.loadState(0);
                runTime.setDeltaT(nextDeltaT);
                
                // First half-step
                Info << "\n---CRAFTS--- Fine steps"
                    << endl;
                if (model.doubleFineStep() != 0)
                {
                    break;
                }

                // Update error scaling factors
                model.calculateScales();

                // Calculate next best time
                oldNextDeltaT = nextDeltaT;
                nextDeltaT = model.calculateNextBestDeltaT();

                // Calculate the residual
                if (model.testConvergence() <= 0)
                {
                    solutionAcceptable = true;
                }
                else
                {
                    nextDeltaT = oldNextDeltaT;
                }

                // Made it through all the steps
                badStep = false;

            } while (false);

            if (badStep || !solutionAcceptable)
            {
                // The solution is unacceptable
                
                if (!model.atsUseAts())
                {
                    FatalErrorIn("craftsPlcFoam::solve")
                        << "Timestep failure.  Cannot continue without "
                        << "adaptive timestepping."
                        << abort(FatalError);
                }
                
                // Now the timestep *must* be reduced by a minimum value
                // i.e. atsMinReductionFactor_ applies
                // First, undouble deltaT, if applicable
                if (runTime.doubled())
                {
                    runTime.halveDeltaT();
                }

                scalar maxDeltaT
                (
                    runTime.deltaT().value()
                  * model.atsMinReductionFactor()
                );
                nextDeltaT = min(nextDeltaT, maxDeltaT);

                // Return to the restore point
                model.loadState(0);

                // Set the new timestep
                scalar oldDeltaT(runTime.deltaT().value());
                runTime.setDeltaTCoarse(nextDeltaT);
                
                // Check for unchanged timestep
                if (runTime.deltaT().value() == oldDeltaT)
                {
                    FatalErrorIn("craftsPlcFoam::solve")
                        << "Required timestep " << nextDeltaT << " is out of "
                        << "limits, " << runTime.deltaT().value()
                        << abort(FatalError);
                }

                // admTime may have limited deltaT, assign the new value
                nextDeltaT = runTime.deltaT().value();
                
                Info << "---CRAFTS--- Reducing timestep to "
                    << runTime.deltaT().value() << " and retrying." << endl;
            }
            else
            {
                scalar oldDeltaT(runTime.deltaT().value());
                
                // Tell the model the good news
                model.transitionToNextTimestep();
                
                runTime.setDeltaTCoarse(nextDeltaT);
                nextDeltaT = runTime.deltaT().value();
                Info << "---CRAFTS--- Adaptive timestepping "
                    << "adjusting next timestep from " << oldDeltaT << " to "
                    << nextDeltaT << endl;
            }
        } // end of adaptive timestep loop
        
        if (runTime.outputTime())
        {
            Info << "---CRAFTS--- Step successful, writing output"
                << endl;
        }
        else
        {
            Info << "---CRAFTS--- Step successful, no output "
                << "generated" << endl;
        }
        model.admVars().updateAverages();
        control.updateAndWrite(runTime);
    }
    
    Info << "---CRAFTS--- End" << endl;
}

// ************************************************************************* //
