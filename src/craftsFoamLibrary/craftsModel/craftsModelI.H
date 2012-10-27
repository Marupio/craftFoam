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

template<int matrixSize>
inline const volVectorField& craftsModel<matrixSize>::U() const
{
    return U_;
}


template<int matrixSize>
inline const surfaceScalarField& craftsModel<matrixSize>::phi() const
{
    return phi_;
}


template<int matrixSize>
inline const word& craftsModel<matrixSize>::functionHooksName() const
{
    return functionHooksName_;
}


template<int matrixSize>
inline autoPtr<admFcHooks<matrixSize> >& craftsModel<matrixSize>::hooks()
{
    return hooks_;
}


template<int matrixSize>
inline const autoPtr<admFcHooks<matrixSize> >&
    craftsModel<matrixSize>::hooks() const
{
    return hooks_;
}


template<int matrixSize>
inline const admFlow& craftsModel<matrixSize>::flow() const
{
    return flow_();
}


template<int matrixSize>
inline admFlow& craftsModel<matrixSize>::flow()
{
    return flow_();
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputOverallSummary() const
{
    return outputOverallSummary_;
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputReactionResidualSummary()
    const
{
    return outputReactionResidualSummary_;
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputReactionResidualDetails()
    const
{
    return outputReactionResidualDetails_;
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputReactionSolverPerformance()
    const
{
    return outputReactionSolverPerformance_;
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputReactionTimestepEstimate()
    const
{
    return outputReactionTimestepEstimate_;
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputAutoSolvePerformance() const
{
    return outputAutoSolvePerformance_;
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputImplicitLoopSummary() const
{
    return outputImplicitLoopSummary_;
}


template<int matrixSize>
inline const bool& craftsModel<matrixSize>::outputImplicitLoopDetails() const
{
    return outputImplicitLoopDetails_;
}


template<int matrixSize>
inline const label& craftsModel<matrixSize>::outerLoopMaxIterations() const
{
    return outerLoopMaxIterations_;
}


template<int matrixSize>
inline const label& craftsModel<matrixSize>::innerLoopMaxIterations() const
{
    return innerLoopMaxIterations_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::standardConvergence() const
{
    return standardConvergence_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::implicitConvergence() const
{
    return implicitConvergence_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::standardScale() const
{
    return standardScale_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::implicitScale() const
{
    return implicitScale_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::standardResidual() const
{
    return standardResidual_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::implicitResidual() const
{
    return implicitResidual_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::implicitDdtResidual() const
{
    return implicitDdtResidual_;
}


template<int matrixSize>
inline const scalarField& craftsModel<matrixSize>::udfDeltaResidual() const
{
    return udfDeltaResidual_;
}


template<int matrixSize>
inline const Switch& craftsModel<matrixSize>::atsUseAts() const
{
    return atsUseAts_;
}


template<int matrixSize>
inline const scalar& craftsModel<matrixSize>::atsConvergenceFactor() const
{
    return atsConvergenceFactor_;
}


template<int matrixSize>
inline const scalar& craftsModel<matrixSize>::atsOverclockFactor() const
{
    return atsOverclockFactor_;
}


template<int matrixSize>
inline const scalar& craftsModel<matrixSize>::atsMaxIncreaseFactor() const
{
    return atsMaxIncreaseFactor_;
}


template<int matrixSize>
inline const scalar& craftsModel<matrixSize>::atsMaxReductionFactor() const
{
    return atsMaxReductionFactor_;
}


template<int matrixSize>
inline const scalar& craftsModel<matrixSize>::atsMinReductionFactor() const
{
    return atsMinReductionFactor_;
}


template<int matrixSize>
inline const Switch& craftsModel<matrixSize>::atsPfUsePf() const
{
    return atsPfUsePf_;
}


template<int matrixSize>
inline const label& craftsModel<matrixSize>::atsPfMeasure() const
{
    return atsPfMeasure_;
}


template<int matrixSize>
inline const scalar& craftsModel<matrixSize>::atsPfBias() const
{
    return atsPfBias_;
}


template<int matrixSize>
inline const label& craftsModel<matrixSize>::lastStepIterations() const
{
    return lastStepIterations_;
}


template<int matrixSize>
inline label& craftsModel<matrixSize>::lastStepIterations()
{
    return lastStepIterations_;
}


template<int matrixSize>
inline const scalar& craftsModel<matrixSize>::lastStepCpuTime() const
{
    return lastStepCpuTime_;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
