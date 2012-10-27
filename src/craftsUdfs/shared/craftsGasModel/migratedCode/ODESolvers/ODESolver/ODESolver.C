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

#include "ODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::ODESolver, 0);
namespace Foam
{
    defineRunTimeSelectionTable(ODESolver, ODE);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODESolver::ODESolver(ODE& ode)
:
    ode_(ode),
    maxStep_(10000),
    yScale_(ode.nEqns()),
    dydx_(ode.nEqns())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::ODESolver::solve
(
    const scalar xStart,
    const scalar xEnd,
    const scalar eps,
    scalar& hEst,
    label& itersDone
) const
{
    scalar x = xStart;
    scalar h = hEst;

    const label nEqns = ode_.nEqns();
    scalarField& y = ode_.coeffs();
    scalarField yOld(y);

    for (label nStep = 0; nStep < maxStep_; nStep++)
    {
        ode_.derivatives(x, y, dydx_);

        for (label i=0; i<nEqns; i++)
        {
            yScale_[i] = mag(y[i]) + mag(dydx_[i]*h) + SMALL;
        }

        if ((x + h - xEnd)*(x + h - xStart) > 0.0)
        {
            h = xEnd - x;
        }

        scalar hNext, hDid;
        { // Scope to delete temporary variables yOld and xOld
            scalarField yOld(y);
            scalar xOld(x);

            solve(x, y, dydx_, eps, yScale_, h, hDid, hNext);
            
            // Trapezoidal rule
            ode_.integralRTSgDt() +=
                (yOld + y) / 2.0 * (x - xOld)
              * ode_.R() * ( 1.0 / 2.0 * (ode_.T(x) + ode_.T(xOld)));
        }

        if ((x - xEnd)*(xEnd - xStart) >= 0.0)
        {
            hEst = hNext;

            // Solution completed.  Update ODE
            ode_.update(xEnd - xStart);
            itersDone = nStep;
            return 0;
        }

        h = hNext;
    }

    scalar percentComplete
    (
        mag((x - xStart) / (xEnd - xStart) * 100)
    );

    WarningIn
    (
        "ODESolver::solve"
        "(const scalar xStart, const scalar xEnd,"
        "scalarField& yStart, const scalar eps, scalar& hEst) const"
    )   << "Too many integration steps after completing " << percentComplete
        << " percent." << endl;
    itersDone = maxStep_;
    return -1;
}


// ************************************************************************* //
