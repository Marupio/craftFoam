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
    Foam::craftsGasFlowOde

Description
    An ODE for solving gas concentrations simultaneously 
    
SourceFiles
    craftsGasFlowOdeI.H
    craftsGasFlowOde.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef craftsGasFlowOde_H
#define craftsGasFlowOde_H

#include "admCoefficientManager.H"
#include "admVariableManager.H"
#include "ODE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

template<int>
class craftsModel;

/*---------------------------------------------------------------------------*\
                        Class craftsGasFlowOde Declaration
\*---------------------------------------------------------------------------*/

template<int matrixSize>
class craftsGasFlowOde
:
    public ODE
{

    //- Reference to model
    craftsModel<matrixSize>& model_;

    // Full geometric fields (these are temporally interpolated)
    
        //- Pointer list to liquid concentrations
        UPtrList<const volScalarField>& slFields_;
        
    // Mesh-related
    
        //- Mesh
        const fvMesh& mesh_;
        
        //- surfacePatch index of gas transfer boundary surface
        const label surfacePatchIndex_;
        
    // Universal constants

        //- Start time
        const scalar tStart_;

        //- End time
        const scalar tEnd_;

        //- Gas volume
        const scalar Vg_;
        
        //- Gas constant
        const scalar R_;

        //- Pipe resistance coefficient
        const scalar kp_;
        
        //- Logistic function sharpness (analytical approximation to step fn)
        const scalar kGasCutOff_;
        
        //- Logistic function offset (gas over-pressure to prevent backflow)
        const scalar pGasOverPressure_;

    // Temporally interpolated variables
    
        //- Gas temperature (at tStart)
        const scalar TgasStart_;

        //- Gas temperature (at tEnd)
        const scalar TgasEnd_;

        //- Atmospheric pressure (at tStart)        
        const scalar pAtmStart_;
        
        //- Atmospheric pressure (at tEnd)
        const scalar pAtmEnd_;
        
        //- Bulk fluid vapour pressure (at tStart)        
        const scalar pLiqStart_;

        //- Bulk fluid vapour pressure (at tEnd)
        const scalar pLiqEnd_;

    // Array variables

        //- Overall diffusion gas transfer coefficient
        const scalarField& kla_;

        //- Henry's coefficient (at tStart)
        const PtrList<scalarField>& khStart_;

        //- Henry's coefficient (at tEnd)
        const PtrList<scalarField>& khEnd_;

        //- Conversion to moles constants
        const scalarField& n_;
        
        //- Bubble source term
        const scalarField& b_;

        //- Upper limit for the coefficients        
        const scalarField& upperLimits_;

        //- Lower limit for the coefficients        
        const scalarField& lowerLimits_;
        
        //- Gas concentration (these are the ode's "coefficients")
        scalarField& sg_;
        
        //- Integration of R* T * sg_ dt - this is just an integration sum that
        //  is used by other objects for accurate mass-transfer calculations
        scalarField integralRTSgDt_;

    // Locally stored variables, derived from the above inputs
    
        //- Average liquid concentrations
        scalarList meanSl_;
        scalarList meanSlOld_;
        
        //- Average Henry's constant
        scalarList meanKh_;
        scalarList meanKhOld_;
        
        //- Sum of the volume of all control volumes at the surface
        scalar surfaceVolume_;
        scalar surfaceVolumeOld_;

public:

    //- Runtime type information
    TypeName("craftsGasFlowOde");


    // Constructors

        //- Construct from components
        craftsGasFlowOde
        (
            craftsModel<matrixSize>& model,
            UPtrList<const volScalarField>& slFields,
            const label surfacePatchIndex,
            const scalar tStart,
            const scalar tEnd,
            const scalar Vg,
            const scalar R,
            const scalar kp,
            const scalar kGasCutOff,
            const scalar pGasOverPressure,
            const scalar TgasStart,
            const scalar TgasEnd,
            const scalar pAtmStart,
            const scalar pAtmEnd,
            const scalar pLiqStart,
            const scalar pLiqEnd,
            const scalarField& kla,
            const PtrList<scalarField>& khStart,
            const PtrList<scalarField>& khEnd,
            const scalarField& n,
            const scalarField& b,
            const scalarField& upperLimits,
            const scalarField& lowerLimits,
            scalarField& sg
        );

    // Destructor

        virtual ~craftsGasFlowOde()
        {}


    // Member functions

        // Access
        
            //- Gas constant
            virtual const scalar& R() const
            {
                return R_;
            }
            
            //- Interpolated temperature
            virtual scalar T(const scalar& currentTime) const;
            
            //- Integral of R * T * sg dt
            virtual scalarField& integralRTSgDt()
            {
                return integralRTSgDt_;
            }
            virtual const scalarField& integralRTSgDt() const
            {
                return integralRTSgDt_;
            }
    
        // Miscellaneous

            //- Perform temporal interpolation with the given start and end
            //  variables
            virtual scalar localInterpolate
            (
                const scalar& start,
                const scalar& end,
                const scalar& currentTime
            ) const;

        // Required by ODE

            //- Return number of equations
            virtual label nEqns() const;

            //- Return access to coefficients
            virtual scalarField& coeffs();

            //- Return reference to coefficients
            virtual const scalarField& coeffs() const;

            //- Apply limits to the coefficients
            virtual void applyCoeffLimits(scalarField& y);

            //- Evaluate derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Evaluate Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta);

    // Friend operators

        // Ostream operator

        template<int mSize>
        friend Ostream& operator<<
        (
            Ostream&,
            const craftsGasFlowOde<mSize>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "craftsGasFlowOde.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
