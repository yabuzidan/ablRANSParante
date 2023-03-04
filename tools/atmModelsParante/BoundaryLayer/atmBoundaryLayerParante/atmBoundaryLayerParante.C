/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "atmBoundaryLayerParante.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::atmBoundaryLayerParante::kappaDefault_ = 0.41;

const Foam::scalar Foam::atmBoundaryLayerParante::CmuDefault_ = 0.09;

const Foam::scalar Foam::atmBoundaryLayerParante::C1Default_ = 0.0;

const Foam::scalar Foam::atmBoundaryLayerParante::C2Default_ = 1.0;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::atmBoundaryLayerParante::init()
{
    if (mag(flowDir_) < small || mag(zDir_) < small)
    {
        FatalErrorInFunction
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    // Ensure direction vectors are normalized
    flowDir_ /= mag(flowDir_);
    zDir_ /= mag(zDir_);

    // (derived from RH:Eq. 6, HW:Eq. 5)
    Ustar_ = kappa_*Uref_/(log((Zref_ + z0_)/z0_));
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmBoundaryLayerParante::atmBoundaryLayerParante()
:
    flowDir_(Zero),
    zDir_(Zero),
    kappa_(0.41),
    Cmu_(0.09),
    C1_(0.0),
    C2_(1.0),
    Uref_(0),
    Zref_(0),
    z0_(0),
    zGround_(0),
    Ustar_(0),
    offset_(false),
    Ulower_(0),
    kLower_(0),
    epsilonLower_(0)
{}


Foam::atmBoundaryLayerParante::atmBoundaryLayerParante
(
    const vector& flowDir,
    const vector& zDir,
    const scalar Uref,
    const scalar Zref,
    const scalarField& z0,
    const scalarField& zGround,
    const scalar kappa,
    const scalar Cmu,
    const scalar C1,
    const scalar C2,
    const scalar Ulower,
    const scalar kLower,
    const scalar epsilonLower
)
:
    flowDir_(flowDir),
    zDir_(zDir),
    kappa_(kappa),
    Cmu_(Cmu),
    C1_(C1),
    C2_(C2),
    Uref_(Uref),
    Zref_(Zref),
    z0_(z0),
    zGround_(zGround),
    Ustar_(z0.size()),
    offset_(Ulower != 0),
    Ulower_(Ulower),
    kLower_(kLower),
    epsilonLower_(epsilonLower)
{
    init();
}


Foam::atmBoundaryLayerParante::atmBoundaryLayerParante
(
    const vectorField& p,
    const dictionary& dict
)
:
    flowDir_(dict.lookup("flowDir")),
    zDir_(dict.lookup("zDir")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", kappaDefault_)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", CmuDefault_)),
    C1_(dict.lookupOrDefault<scalar>("C1", C1Default_)),
    C2_(dict.lookupOrDefault<scalar>("C2", C2Default_)),
    Uref_(dict.lookup<scalar>("Uref")),
    Zref_(dict.lookup<scalar>("Zref")),
    z0_("z0", dict, p.size()),
    zGround_("zGround", dict, p.size()),
    Ustar_(p.size()),
    offset_(dict.found("Ulower")),
    Ulower_(dict.lookupOrDefault<scalar>("Ulower", 0)),
    kLower_(dict.lookupOrDefault<scalar>("kLower", 0)),
    epsilonLower_(dict.lookupOrDefault<scalar>("epsilonLower", 0))
{
    init();
}


Foam::atmBoundaryLayerParante::atmBoundaryLayerParante
(
    const atmBoundaryLayerParante& abl,
    const fvPatchFieldMapper& mapper
)
:
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    C1_(abl.C1_),
    C2_(abl.C2_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    z0_(mapper(abl.z0_)),
    zGround_(mapper(abl.zGround_)),
    Ustar_(mapper(abl.Ustar_)),
    offset_(abl.offset_),
    Ulower_(abl.Ulower_),
    kLower_(abl.kLower_),
    epsilonLower_(abl.epsilonLower_)
{}


Foam::atmBoundaryLayerParante::atmBoundaryLayerParante(const atmBoundaryLayerParante& abl)
:
    flowDir_(abl.flowDir_),
    zDir_(abl.zDir_),
    kappa_(abl.kappa_),
    Cmu_(abl.Cmu_),
    C1_(abl.C1_),
    C2_(abl.C2_),
    Uref_(abl.Uref_),
    Zref_(abl.Zref_),
    z0_(abl.z0_),
    zGround_(abl.zGround_),
    Ustar_(abl.Ustar_),
    offset_(abl.offset_),
    Ulower_(abl.Ulower_),
    kLower_(abl.kLower_),
    epsilonLower_(abl.epsilonLower_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmBoundaryLayerParante::autoMap(const fvPatchFieldMapper& m)
{
    m(z0_, z0_);
    m(zGround_, zGround_);
    m(Ustar_, Ustar_);
}


void Foam::atmBoundaryLayerParante::rmap
(
    const atmBoundaryLayerParante& blptf,
    const labelList& addr
)
{
    z0_.rmap(blptf.z0_, addr);
    zGround_.rmap(blptf.zGround_, addr);
    Ustar_.rmap(blptf.Ustar_, addr);
}


Foam::tmp<Foam::vectorField> Foam::atmBoundaryLayerParante::U
(
    const vectorField& p
) const
{
    const scalarField Un
    (
        // (Parente:Eq 5., YGCJ:Table 1, RH:Eq. 6, HW:Eq. 5)
        (Ustar_/kappa_)*log(max((zDir_ & p) - zGround_ + z0_, z0_)/z0_)
    );

    if (offset_)
    {
        return flowDir_*Un + flowDir_*Ulower_;
    }
    else
    {
        return flowDir_*Un;
    }
}


Foam::tmp<Foam::scalarField> Foam::atmBoundaryLayerParante::k
(
    const vectorField& p
) const
{
    tmp<scalarField> tk
    (
        // (Parente:Eq. 15)
        C1_*log(max((zDir_ & p) - zGround_ + z0_, z0_)) + C2_
        
    //     // (YGCJ:Eq. 21; RH:Eq. 7, HW:Eq. 6 when C1=0 and C2=1)
    //     sqr(Ustar_)/sqrt(Cmu_)
    //    *sqrt(C1_*log(max((zDir_ & p) - zGround_ + z0_, z0_)/z0_) + C2_)
    );

    if (offset_)
    {
        const scalarField z((zDir_ & p) - zGround_);
        tk.ref() = pos0(z)*tk() + neg(z)*kLower_;
    }

    return tk;
}


Foam::tmp<Foam::scalarField> Foam::atmBoundaryLayerParante::epsilon
(
    const vectorField& p
) const
{
    tmp<scalarField> tepsilon
    (
        // (Parente:Eq. 7; RH:Eq. 8, HW:Eq. 7)
        pow3(Ustar_)/(kappa_*((zDir_ & p) - zGround_ + z0_))

    //     // (YGCJ:Eq. 22; RH:Eq. 8, HW:Eq. 7 when C1=0 and C2=1)
    //     pow3(Ustar_)/(kappa_*((zDir_ & p) - zGround_ + z0_))
    //    *sqrt(C1_*log(max((zDir_ & p) - zGround_ + z0_, z0_)/z0_) + C2_)
    );

    if (offset_)
    {
        const scalarField z((zDir_ & p) - zGround_);
        tepsilon.ref() = pos0(z)*tepsilon() + neg(z)*epsilonLower_;
    }

    return tepsilon;
}


void Foam::atmBoundaryLayerParante::write(Ostream& os) const
{
    writeEntry(os, "z0", z0_) ;
    writeEntry(os, "flowDir", flowDir_);
    writeEntry(os, "zDir", zDir_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "Cmu", Cmu_);
    writeEntry(os, "C1", C1_);
    writeEntry(os, "C2", C2_);
    writeEntry(os, "Uref", Uref_);
    writeEntry(os, "Zref", Zref_);

    if (offset_)
    {
        writeEntry(os, "Ulower", Ulower_);
        writeEntry(os, "kLower", kLower_);
        writeEntry(os, "epsilonLower", epsilonLower_);
    }

    writeEntry(os, "zGround", zGround_);
}


// ************************************************************************* //
