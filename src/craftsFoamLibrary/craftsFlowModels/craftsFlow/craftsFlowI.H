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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const admReactionReader& craftsFlow::model() const
{
    return model_;
}


inline const admTime& craftsFlow::runTime() const
{
    return runTime_;
}


inline const fvMesh& craftsFlow::mesh() const
{
    return mesh_;
}


inline const volVectorField& craftsFlow::U() const
{
    return UGlobal_;
}


inline volVectorField& craftsFlow::U()
{
    return UGlobal_;
}


inline const surfaceScalarField& craftsFlow::phi() const
{
    return phiGlobal_;
}


inline surfaceScalarField& craftsFlow::phi()
{
    return phiGlobal_;
}


inline const volScalarField& craftsFlow::p() const
{
    return pGlobal_;
}


inline volScalarField& craftsFlow::p()
{
    return pGlobal_;
}


inline bool craftsFlow::atsEnableWithFlowModel() const
{
    return atsEnableWithFlowModel_;
}


inline const wordList& craftsFlow::atsIgnore() const
{
    return atsIgnore_;
}


inline bool craftsFlow::atsIgnorePhi() const
{
    return atsIgnorePhi_;
}


inline bool craftsFlow::atsIgnoreP() const
{
    return atsIgnoreP_;
}


inline const scalar& craftsFlow::atsPhiTolerance() const
{
    return atsPhiTolerance_;
}


inline const scalar& craftsFlow::atsPTolerance() const
{
    return atsPTolerance_;
}


inline const scalar& craftsFlow::atsPhiMinScale() const
{
    return atsPhiMinScale_;
}


inline const scalar& craftsFlow::atsPMinScale() const
{
    return atsPMinScale_;
}


inline const scalar& craftsFlow::atsPhiScale() const
{
    return atsPhiScale_;
}


inline const scalar& craftsFlow::atsPScale() const
{
    return atsPScale_;
}


inline const bool& craftsFlow::subStepping() const
{
    return subStepping_;
}


inline const craftsFlow::enumSubSteppingType&
    craftsFlow::subSteppingType() const
{
    return subSteppingType_;
}


inline const label& craftsFlow::nSubSteps() const
{
    return ssCurrentNSubSteps_;
}


inline const bool& craftsFlow::outputFlowSolverPerformance() const
{
    return outputFlowSolverPerformance_;
}


inline const bool& craftsFlow::outputFlowResiduals() const
{
    return outputFlowResiduals_;
}


inline const bool& craftsFlow::outputFlowTimestepEstimate() const
{
    return outputFlowTimestepEstimate_;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
