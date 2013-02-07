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

#include "craftsModel.H"
#include "IOstreams.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::implicitAutoSolve
(
    admImplicitVariable& impVar,
    label cellIndex
)
{
    // Same as implicitAutoSolveField, except for only one control volume.
    // See comments at the start of that function.  This one performs full
    // field fvMatrix calculations, and throws away all cells except for
    // cellIndex.  Inefficient - implicitAutoSolveField is preferred.
    if (!impVar.autoSolve())
    {
        FatalErrorIn("craftsModel<matrixSize>::implicitAutoSolve")
            << "implicitAutoSolve called for a variable that does not "
            << "have autoSolve set to 'yes' in the admVariableDict."
            << abort(FatalError);
    }

    scalar& Svar(impVar().internalField()[cellIndex]);

    label maxIter(impVar.autoSolveMaxIter());
    scalar convergence(impVar.autoSolveConvergence());
    label autoSolveIndex(impVar.autoSolveIndex());
    label varIndex(impVar.localIndex());

    label iterations(0);

    for (; iterations < maxIter; iterations++)
    {
        scalar Tb(0.0);
        scalar Td(0.0);
        scalar sumTiSi(0.0);
        scalar R(0.0);
        scalar dRdS(0.0);
        const scalar Sold(impVar.evaluate(cellIndex));

        // Calculate Tb and Td - create transport matrix
        fvScalarMatrix SEqn
        (
            fvm::ddt(impVar())
          + fvm::div
            (
                phi_, impVar()
            )
          - fvm::laplacian
            (
                impVar.gamma(),
                impVar()
            )
        );
        
        // Absorb boundary conditions (stolen from fvMatrix<scalar>::solve)
        SEqn.addBoundaryDiag(SEqn.diag(), 0);
        scalarField totalSource = SEqn.source();
        SEqn.addBoundarySource(totalSource, false);
        
        // Add Td - diagonal transport term
        Td = SEqn.diag()[cellIndex];

        // Calculate Tb - transport source term contribution
        Tb += SEqn.source()[cellIndex];

        // Calculate Tb - transport owners contribution
        forAll(transportOwners_[cellIndex], i)
        {
            label faceI(transportOwners_[cellIndex][i]);
            label cellI(transportOwnersCells_[cellIndex][i]);
            sumTiSi += 
                SEqn.lower()[faceI] * impVar.evaluate(cellI);
        }
        
        // Calculate Tb - transport neighbours contribution
        forAll(transportNeighbours_[cellIndex], i)
        {
            label faceI(transportNeighbours_[cellIndex][i]);
            label cellI(transportNeighboursCells_[cellIndex][i]);
            sumTiSi +=
                SEqn.upper()[faceI] * impVar.evaluate(cellI);
        }

        // Calculate dRdS - add termI contributions
        forAll(dRdSAutoSolveTermIk_[autoSolveIndex], termI)
        {
            const admReaction& reaction
            (
                dRdSAutoSolveTermIk_[autoSolveIndex][termI]
            );
            dRdS += reaction.yield(impVar).evaluate(cellIndex)
              * reaction.rate().ddy(impVar, cellIndex);
        }

        // Calculate dRdS - add term II contributions
        forAll(dRdSAutoSolveTermIIk_[autoSolveIndex], termII)
        {
            const admReaction& reaction
            (
                dRdSAutoSolveTermIIk_[autoSolveIndex][termII]
            );
            dRdS += reaction.yield(impVar).ddy(impVar, cellIndex)
              * reaction.rate().evaluate(cellIndex);
        }

        // Calculate R
        forAll(nonZeroImplicitYields_[varIndex], yieldIndex)
        {
            const admReaction& reaction
            (
                nonZeroImplicitYields_[varIndex][yieldIndex]
            );

            R += reaction.yield(impVar).evaluate(cellIndex)
          * reaction.rate().evaluate(cellIndex);
        }

        // Calculate next value for S
        Svar = (Tb + R - dRdS * Svar - sumTiSi) / (Td - dRdS);

        // Test convergence
        if
        (
            mag
            (
                (Svar - Sold) / implicitScale_[impVar.localIndex()]
            )  < convergence
        )
        {
            break;
        }
    }

    // Apply restrictions
    admVars_.applyLimits(impVar, cellIndex);
    
    if (iterations == maxIter)
    {
        WarningIn("craftsModel<matrixSize>::implicitAutoSolve")
            << "Exceeded max iterations for variable " << impVar.name()
            << " at cell index " << cellIndex
            << endl;

        // Add iterations to performance stats
        lastStepIterations_ += iterations;

        return -1;
    }
    
    // Add iterations to performance stats
    lastStepIterations_ += iterations;
    
    return 0;
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::implicitAutoSolveField
(
    admImplicitVariable& impVar
)
{
    // Solves AS = B, where
    //  A is the diagonal coefficient, given by:
    //      A = Td - dR/dS,
    //      where:
    //          Td is the diagonal term from the transport equation, and
    //          R is the reaction source term.
    //  S is the variable, and
    //  B is the matrix source term, given by:
    //      B = Tb - sum(Ti Si) + R - dR/dS S*,
    //      where:
    //          Tb is the matrix source term from the transport equation,
    //          Ti are the owner/neighbour terms from the transport equation,
    //          Si are the S value in the owner/neighbour cell, and
    //          S* is the previous value of S.
    //  In other words:
    //      S = (Tb + R - dR/dS S* - sum(Ti Si)) / (Td - dR/dS)
    if (!impVar.autoSolve())
    {
        FatalErrorIn("craftsModel<matrixSize>::implicitAutoSolveField")
            << "implicitAutoSolve called for a variable that does not "
            << "have autoSolve set to 'yes' in the admVariableDict."
            << abort(FatalError);
    }

    scalarField& Svar(impVar().internalField());

    label maxIter(impVar.autoSolveMaxIter());
    scalar convergence(impVar.autoSolveConvergence());
    scalar convergedTo(VGREAT);
    label autoSolveIndex(impVar.autoSolveIndex());
    label varIndex(impVar.localIndex());

    label iterations(0);

    for (; iterations < maxIter; iterations++)
    {
        scalarField Tb(mesh_.nCells(), 0.0);
        scalarField Td(mesh_.nCells(), 0.0);
        scalarField sumTiSi(mesh_.nCells(), 0.0);
        scalarField R(mesh_.nCells(), 0.0);
        scalarField dRdS(mesh_.nCells(), 0.0);
        const scalarField& Sold(impVar.evaluateField());

        // Calculate Tb and Td - create transport matrix
        fvScalarMatrix SEqn
        (
            fvm::ddt(impVar())
          + fvm::div
            (
                phi_, impVar()
            )
          - fvm::laplacian
            (
                impVar.gamma(),
                impVar()
            )
        );
        
        // Absorb boundary conditions (stolen from fvMatrix<scalar>::solve)
        SEqn.addBoundaryDiag(SEqn.diag(), 0);
        scalarField totalSource = SEqn.source();
        SEqn.addBoundarySource(totalSource, false);
        
        // Add Td - diagonal transport term
        Td = SEqn.diag();

        // Calculate Tb - transport source term contribution
        Tb += SEqn.source();

        // Calculate Tb - transport owners contribution
        // This is a hack job and can be improved
        forAll(transportOwners_, cellIndex)
        {
            forAll(transportOwners_[cellIndex], i)
            {
                label faceI(transportOwners_[cellIndex][i]);
                label cellI(transportOwnersCells_[cellIndex][i]);
                sumTiSi[cellIndex] += 
                    SEqn.lower()[faceI] * impVar.evaluate(cellI);
            }
        }
        
        // Calculate Tb - transport neighbours contribution
        // This is a hack job and can be improved
        forAll(transportNeighbours_, cellIndex)
        {
            forAll(transportNeighbours_[cellIndex], i)
            {
                label faceI(transportNeighbours_[cellIndex][i]);
                label cellI(transportNeighboursCells_[cellIndex][i]);
                sumTiSi[cellIndex] +=
                    SEqn.upper()[faceI] * impVar.evaluate(cellI);
            }
        }

        // Calculate dRdS - add termI contributions
        forAll(dRdSAutoSolveTermIk_[autoSolveIndex], termI)
        {
            const admReaction& reaction
            (
                dRdSAutoSolveTermIk_[autoSolveIndex][termI]
            );
            dRdS += reaction.yield(impVar).evaluateField()
              * reaction.rate().ddyField(impVar);
        }

        // Calculate dRdS - add term II contributions
        forAll(dRdSAutoSolveTermIIk_[autoSolveIndex], termII)
        {
            const admReaction& reaction
            (
                dRdSAutoSolveTermIIk_[autoSolveIndex][termII]
            );
            dRdS += reaction.yield(impVar).ddyField(impVar)
              * reaction.rate().evaluateField();
        }

        // Calculate R
        forAll(nonZeroImplicitYields_[varIndex], yieldIndex)
        {
            const admReaction& reaction
            (
                nonZeroImplicitYields_[varIndex][yieldIndex]
            );

            R += reaction.yield(impVar).evaluateField()
          * reaction.rate().evaluateField();
        }

        // Calculate next value for S
        Svar = (Tb + R - dRdS * Svar - sumTiSi) / (Td - dRdS);

        // Test convergence
        convergedTo = calculateRmsError
        (
            Svar,
            Sold,
            implicitScale_[impVar.localIndex()]
        ) - convergence;

        if (convergedTo < 0)
        {
            break;
        }
    }

    // Apply restrictions
    admVars_.applyLimits(impVar);
    
    // Add iterations to performance stats
    lastStepIterations_ += iterations + 1;

    if (outputAutoSolvePerformance_)
    {
        Info << "implicitAutoSolvePerformance: n = " << label(iterations + 1)
            << ", " << impVar.name() << " res = " << convergedTo << endl;
    }

    if (iterations == maxIter)
    {
        WarningIn("craftsModel<matrixSize>::implicitAutoSolveField")
            << "Exceeded max iterations for variable " << impVar.name()
            << endl;

        return -1;
    }
    return 0;
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::implicitAutoSolveAll
(
    label cellIndex
)
{
    label returnMe(0);
    forAll(admVars_.implicitAutoSolve(), impIndex)
    {
        admImplicitVariable& impVar(admVars_.implicitAutoSolve(impIndex));
        label errorLevel(implicitAutoSolve(impVar, cellIndex));
        switch (errorLevel)
        {
            default:
                // Unknown return value, assume total failure
            case -2:
                // Total failure, stop immediately
                return errorLevel;
            case -1:
                returnMe = errorLevel;
                break;
            case 0:
                // All is okay
                break;
        }
    }
    return returnMe;
}


template <int matrixSize>
Foam::label Foam::craftsModel<matrixSize>::implicitAutoSolveAllField()
{
    label returnMe(0);
    forAll(admVars_.implicitAutoSolve(), impIndex)
    {
        admImplicitVariable& impVar(admVars_.implicitAutoSolve(impIndex));
        label errorLevel(implicitAutoSolveField(impVar));
        switch (errorLevel)
        {
            default:
                // Unknown return value, assume total failure
            case -2:
                // Total failure, stop immediately
                return errorLevel;
            case -1:
                returnMe = errorLevel;
                break;
            case 0:
                // All is okay
                break;
        }
    }
    return returnMe;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
