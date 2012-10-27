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

Class
    Foam::ODE

Description
    Abstract base class for the ODE solvers.

\*---------------------------------------------------------------------------*/

#ifndef ODE_H
#define ODE_H

#include "scalarField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ODE Declaration
\*---------------------------------------------------------------------------*/

class ODE
{
public:

    // Constructors

        //- Construct null
        ODE()
        {}


    // Destructor

        virtual ~ODE()
        {}


    // Member Functions

        // Hack job functions
        
            //- Gas constant
            virtual const scalar& R() const = 0;
            
            //- Interpolated temperature
            virtual scalar T(const scalar& currentTime) const = 0;
            
            //- Integral of R * T * sg dt
            virtual scalarField& integralRTSgDt() = 0;
            virtual const scalarField& integralRTSgDt() const = 0;
    
            //- Perform temporal interpolation with the given start and end
            //  variables
            virtual scalar localInterpolate
            (
                const scalar& start,
                const scalar& end,
                const scalar& currentTime
            ) const = 0;

        // Standard ODE functions

            //- Return number of equations
            virtual label nEqns() const = 0;

            //- Return access to coefficients
            virtual scalarField& coeffs() = 0;

            //- Return reference to coefficients
            virtual const scalarField& coeffs() const = 0;

            //- Apply limits to the coefficients
            virtual void applyCoeffLimits(scalarField& coeffs)
            {
                // Default behaviour is to do nothing as this is fastest
            }

            //- Evaluate derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const = 0;

            //- Evaluate Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const = 0;

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
