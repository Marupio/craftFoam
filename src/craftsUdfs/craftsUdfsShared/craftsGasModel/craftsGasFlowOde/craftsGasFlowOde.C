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

#include "craftsGasFlowOde.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::craftsGasFlowOde<matrixSize>::craftsGasFlowOde
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
)
:
    model_(model),
    slFields_(slFields),
    mesh_(slFields_[0].mesh()),
    surfacePatchIndex_(surfacePatchIndex),
    tStart_(tStart),
    tEnd_(tEnd),
    Vg_(Vg),
    R_(R),
    kp_(kp),
    kGasCutOff_(kGasCutOff),
    pGasOverPressure_(pGasOverPressure),
    TgasStart_(TgasStart),
    TgasEnd_(TgasEnd),
    pAtmStart_(pAtmStart),
    pAtmEnd_(pAtmEnd),
    pLiqStart_(pLiqStart),
    pLiqEnd_(pLiqEnd),
    kla_(kla),
    khStart_(khStart),
    khEnd_(khEnd),
    n_(n),
    b_(b),
    upperLimits_(upperLimits),
    lowerLimits_(lowerLimits),
    sg_(sg),
    integralRTSgDt_(sg.size(), 0.0),
    meanSl_(slFields.size()),
    meanSlOld_(slFields.size()),
    meanKh_(slFields.size()),
    meanKhOld_(slFields.size()),
    surfaceVolume_(0.0),
    surfaceVolumeOld_(0.0)
{
    // Calculate locally-stored derived values
    forAll(meanSl_, i)
    {
        meanSl_[i] = slFields_[i].weightedAverage(mesh_.V()).value();
        meanSlOld_[i] =
            slFields_[i].oldTime().weightedAverage(mesh_.V()).value();
        meanKh_[i] = gSum(mesh_.V() * khEnd_[i]) / gSum(mesh_.V());
        meanKhOld_[i] = gSum(mesh_.V() * khStart_[i]) / gSum(mesh_.V());
    }

    const labelList& surfaceNeighbourIndices
    (
        slFields_[0].boundaryField()[surfacePatchIndex_]
            .patch().faceCells()
    );
    forAll(surfaceNeighbourIndices, listI)
    {
        label cellIndex(surfaceNeighbourIndices[listI]);
        surfaceVolume_ += mesh_.V().field()[cellIndex];
    }

    const labelList& oldSurfaceNeighbourIndices
    (
        slFields_[0].oldTime().boundaryField()[surfacePatchIndex_]
            .patch().faceCells()
    );
    forAll(oldSurfaceNeighbourIndices, listI)
    {
        label cellIndex(oldSurfaceNeighbourIndices[listI]);
        surfaceVolumeOld_ += mesh_.V().field()[cellIndex];
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<int matrixSize>
Foam::scalar Foam::craftsGasFlowOde<matrixSize>::T(const scalar& currentTime) const
{
    return localInterpolate(TgasStart_, TgasEnd_, currentTime);
}


template<int matrixSize>
Foam::scalar Foam::craftsGasFlowOde<matrixSize>::localInterpolate
(
    const scalar& start,
    const scalar& end,
    const scalar& currentTime
) const
{
    return (end - start) / (tEnd_ - tStart_) * (currentTime - tStart_) + start;
}


template<int matrixSize>
Foam::label Foam::craftsGasFlowOde<matrixSize>::nEqns() const
{
    return sg_.size();
}


template<int matrixSize>
Foam::scalarField& Foam::craftsGasFlowOde<matrixSize>::coeffs()
{
    return sg_;
}


template<int matrixSize>
const Foam::scalarField& Foam::craftsGasFlowOde<matrixSize>::coeffs() const
{
    return sg_;
}


template<int matrixSize>
void Foam::craftsGasFlowOde<matrixSize>::applyCoeffLimits(scalarField& y)
{
    forAll(y, i)
    {
        y[i] = max(y[i], lowerLimits_[i]);
        y[i] = min(y[i], upperLimits_[i]);
    }
}


template<int matrixSize>
void Foam::craftsGasFlowOde<matrixSize>::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    scalar Tgas(localInterpolate(TgasStart_, TgasEnd_, x));
    scalar pAtm(localInterpolate(pAtmStart_, pAtmEnd_, x));
    scalar pLiq(localInterpolate(pLiqStart_, pLiqEnd_, x));
    
    forAll(y, i)
    {
        scalar sli
        (
            localInterpolate
            (
                meanSlOld_[i],
                meanSl_[i],
                x
            )
        );
        scalar khi
        (
            localInterpolate
            (
                meanKhOld_[i],
                meanKh_[i],
                x
            )
        );
        scalar Vl
        (
            localInterpolate
            (
                surfaceVolumeOld_,
                surfaceVolume_,
                x
            )
        );

        scalar diffusion
        (
            Vl * (sli - R_ * Tgas * khi * y[i]) * kla_[i] / Vg_
        );

        scalar pGas(0.0);
        forAll(y, j)
        {
            pGas += y[j] / n_[j];
        }
        pGas *= R_ * Tgas;
        pGas += pLiq;
        scalar deltaP(pGas - pAtm);
        scalar Hexp
        (
            -2 * kGasCutOff_ * (deltaP - pGasOverPressure_)
        );

        if (mag(Hexp) >= log(VGREAT / 2))
        {
            Hexp = pos(Hexp) * VGREAT / 2;
        }
        else
        {
            Hexp = exp(Hexp);
        }
        scalar qGas
        (
            kp_ * deltaP * pGas / pAtm / (1 + Hexp)
        );
        dydx[i] = b_[i] + diffusion - y[i] * qGas / Vg_;
    }
}


template<int matrixSize>
void Foam::craftsGasFlowOde<matrixSize>::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdx,
    scalarSquareMatrix& dfdy
) const
{
    scalar Tgas(localInterpolate(TgasStart_, TgasEnd_, x));
    scalar pAtm(localInterpolate(pAtmStart_, pAtmEnd_, x));
    scalar pLiq(localInterpolate(pLiqStart_, pLiqEnd_, x));
    
    scalar pGas(sum(y / n_ * R_ * Tgas) + pLiq);
    scalar deltaP(pGas - pAtm);
    scalar Hexp
    (
        -2 * kGasCutOff_ * (deltaP - pGasOverPressure_)
    );
    if (mag(Hexp) >= log(VGREAT / 2)) 
    {
        Hexp = pos(Hexp) * VGREAT / 2;
    }
    else
    {
        Hexp = exp(Hexp);
    }

    scalar qGas(kp_ * deltaP * pGas / pAtm / (1 + Hexp));

    scalarField dPdy(R_ * Tgas / n_);
    scalarField dHdy
    (
        Hexp / (1 + Hexp) / (1 + Hexp) * 2 * kGasCutOff_ * dPdy
    );
    scalarField dQdy
    (
        kp_ / pAtm * (2 * pGas - pAtm) * dPdy / (1 + Hexp)
      + kp_ * pGas / pAtm * deltaP * dHdy
    );
    
    forAll(y, i)
    {
        scalar khi
        (
            localInterpolate
            (
                meanKhOld_[i],
                meanKh_[i],
                x
            )
        );
        scalar Vl
        (
            localInterpolate
            (
                surfaceVolumeOld_,
                surfaceVolume_,
                x
            )
        );

        forAll(y, j)
        {
            dfdy[i][j] = -y[i] * dQdy[j];
        }
        dfdy[i][i] -= qGas / Vg_ - kla_[i] / Vg_ * Vl * khi * R_ * Tgas;
    }
}


template<int matrixSize>
void Foam::craftsGasFlowOde<matrixSize>::update(const scalar delta)
{}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<int matrixSize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os, const craftsGasFlowOde<matrixSize>& gfo
)
{
    os  << "/*tStart    */" << token::TAB << gfo.tStart_ << nl
        << "/*tEnd      */" << token::TAB << gfo.tEnd_ << nl
        << "/*Vg        */" << token::TAB << gfo.Vg_ << nl
        << "/*R         */" << token::TAB << gfo.R_ << nl
        << "/*kp        */" << token::TAB << gfo.kp_ << nl
        << "/*TgasStart */" << token::TAB << gfo.TgasStart_ << nl
        << "/*TgasEnd   */" << token::TAB << gfo.TgasEnd_ << nl
        << "/*pAtmStart */" << token::TAB << gfo.pAtmStart_ << nl
        << "/*pAtmEnd   */" << token::TAB << gfo.pAtmEnd_ << nl
        << "/*pLiqStart */" << token::TAB << gfo.pLiqStart_ << nl
        << "/*pLiqEnd   */" << token::TAB << gfo.pLiqEnd_ << nl
        << "/*kla       */" << token::TAB << gfo.kla_ << nl
        << "/*khStart   */" << token::TAB << gfo.khStart_ << nl
        << "/*khEnd     */" << token::TAB << gfo.khEnd_ << nl
        << "/*n         */" << token::TAB << gfo.n_ << nl
        << "/*b         */" << token::TAB << gfo.b_ << nl
        << "/*s         */" << token::TAB << gfo.sg_ << endl;

    os.check("Ostream& operator<<(Ostream&, const craftsGasFlowOde&)");

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
