/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::atmNutWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent kinematic viscosity for
    atmospheric velocity profiles.  It is designed to be used in conjunction
    with the atmBoundaryLayerInletVelocity boundary condition.  The values
    are calculated using:

        \f[
            U = frac{U_f}{\kappa} ln(\frac{z + z_0}{z_0})
        \f]

    where

    \vartable
        U_f    | frictional velocity
        \kappa | Von Karman's constant
        z_0    | surface roughness length
        z      | vertical co-ordinate
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        z0           | surface roughness length| yes         |
    \endtable

    References:
    \verbatim
        Theoretical expressions (tags:RH, SBJM, SM):
            Richards, P. J., & Hoxey, R. P. (1993).
            Appropriate boundary conditions for computational wind
            engineering models using the k-ε turbulence model.
            In Computational Wind Engineering 1 (pp. 145-153).
            DOI:10.1016/B978-0-444-81688-7.50018-8

            Sørensen, N. N., Bechmann, A., Johansen, J., Myllerup, L.,
            Botha, P., Vinther, S., & Nielsen, B. S. (2007).
            Identification of severe wind conditions using
            a Reynolds Averaged Navier-Stokes solver.
            In Journal of Physics: Conference
            series (Vol. 75, No. 1, p. 012053).
            DOI:10.1088/1742-6596/75/1/012053

            Sumner, J., & Masson, C. (2012).
            k−ε simulations of the neutral atmospheric boundary layer:
            analysis and correction of discretization errors on practical grids.
            International journal for numerical
            methods in fluids, 70(6), 724-741.
            DOI:10.1002/fld.2709
    \endverbatim

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            atmNutWallFunction;
        z0              uniform 0;
    }
    \endverbatim

See also
    - Foam::nutkWallFunctionFvPatchScalarField
    - Foam::atmNutkWallFunctionFvPatchScalarField
    - Foam::atmNutUWallFunctionFvPatchScalarField

SourceFiles
    atmNutWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef atmNutWallFunctionFvPatchScalarField_H
#define atmNutWallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class atmNutWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmNutWallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- Surface roughness length
        scalarField z0_;


    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("atmNutWallFunction");


    // Constructors

        //- Construct from patch and internal field
        atmNutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmNutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- atmNutWallFunctionFvPatchScalarField
        //- onto a new patch
        atmNutWallFunctionFvPatchScalarField
        (
            const atmNutWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmNutWallFunctionFvPatchScalarField
        (
            const atmNutWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new atmNutWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        atmNutWallFunctionFvPatchScalarField
        (
            const atmNutWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new atmNutWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access functions

            // Return z0
            scalarField& z0()
            {
                return z0_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
