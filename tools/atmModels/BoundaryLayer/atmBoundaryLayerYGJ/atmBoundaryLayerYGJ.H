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

Class
    Foam::atmBoundaryLayerYGJ

Description
    Base class to set log-law type ground-normal inlet boundary conditions for
    wind velocity and turbulence quantities for homogeneous, two-dimensional,
    dry-air, equilibrium and neutral atmospheric boundary layer (ABL) modelling.

    The ground-normal profile expressions are due to \c YGCJ
    (refer to references below) whereat \c RH expressions were generalised:

    \f[
        u = \frac{U^*}{\kappa} \ln \left( \frac{z - d + z_0}{z_0} \right)
    \f]

    
    \f[
        k = \frac{(U^*)^2}{\sqrt{C_mu}}
        \sqrt{C_1 \ln \left( \frac{z - d + z_0}{z_0} \right) + C_2}
    \f]

    \f[
        \epsilon = \frac{(U^*)^3}{\kappa (z - d + z_0)}
        \sqrt{C_1 \ln \left( \frac{z - d + z_0}{z_0} \right) + C_2}
    \f]
    
    \f[
        U^* =
            \frac{u_{ref} \kappa}{\ln\left(\frac{Z_{ref} + z_0}{z_0}\right)}
    \f]

    where
    \vartable
        U^*     | Friction velocity
        \kappa  | von Karman's constant
        C_mu    | Turbulence viscosity coefficient
        z       | Vertical coordinate
        z_0     | Surface roughness height [m]
        z_g     | Minimum z-coordinate [m]
    \endvartable
    and
    
    
    where
    \vartable
        U_{ref} | Reference velocity at \f$Z_{ref}\f$ [m/s]
        Z_{ref} | Reference height [m]
 \endvartable

    Reference:
    \verbatim
        The ground-normal profile expressions (tag:RH):
            Richards, P. J., & Hoxey, R. P. (1993).
            Appropriate boundary conditions for computational wind
            engineering models using the k-ε turbulence model.
            In Computational Wind Engineering 1 (pp. 145-153).
            DOI:10.1016/B978-0-444-81688-7.50018-8

        Modifications to preserve the profiles downstream (tag:HW):
            Hargreaves, D. M., & Wright, N. G. (2007).
            On the use of the k–ε model in commercial CFD software
            to model the neutral atmospheric boundary layer.
            Journal of wind engineering and
            industrial aerodynamics, 95(5), 355-369.
            DOI:10.1016/j.jweia.2006.08.002

        Expression generalisations to allow height
        variation for turbulence quantities (tag:YGCJ):
            Yang, Y., Gu, M., Chen, S., & Jin, X. (2009).
            New inflow boundary conditions for modelling the neutral equilibrium
            atmospheric boundary layer in computational wind engineering.
            J. of Wind Engineering and Industrial Aerodynamics, 97(2), 88-95.
            DOI:10.1016/j.jweia.2008.12.001

        The generalised ground-normal profile expression for omega (tag:YGJ):
            Yang, Y., Gu, M., & Jin, X., (2009).
            New inflow boundary conditions for modelling the
            neutral equilibrium atmospheric boundary layer in SST k-ω model.
            In: The Seventh Asia-Pacific Conference on Wind Engineering,
            November 8-12, Taipei, Taiwan.

        Theoretical remarks (tag:E):
            Emeis, S. (2013).
            Wind Energy Meteorology: Atmospheric
            Physics for Wind Power Generation.
            Springer-Verlag Berlin Heidelberg.
            DOI:10.1007/978-3-642-30523-8
    \endverbatim

Usage
    \table
        Property     | Description                      | Required  | Default
        flowDir      | Flow direction                   | yes       |
        zDir         | Vertical direction               | yes       |
        kappa        | von Karman's constant            | no        | 0.41
        Cmu          | Turbulence viscosity coefficient | no        | 0.09
        Uref         | Reference velocity [m/s]         | yes       |
        Zref         | Reference height [m]             | yes       |
        z0           | Surface roughness height [m]     | yes       |
        zGround      | Minimum z-coordinate [m]         | yes       |
        Ulower       | Velocity below the BL            | no        |
        kLower       | Turbulence k below the BL        | no        |
        epsilonLower | Turbulence epsilon below the BL  | no        |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    ground
    {
        type            atmBoundaryLayerYGJInletVelocity;
        flowDir         (1 0 0);
        zDir            (0 0 1);
        Uref            10.0;
        Zref            20.0;
        z0              uniform 0.1;
        zGround         uniform 0.0;
    }
    \endverbatim

Note
    - The \c RH expressions are special cases of those in \c YGCJ when \c C1=0
    and \c C2=1. Both \c C1 and \c C2 can be determined by nonlinear fitting
    of (\c YGCJ:Eqs. 19-20) with an experimental dataset for \c k. By default,
    \c atmBoundaryLayerInlet boundary conditions compute \c RH expressions.
    - \c z is the ground-normal height relative to the global minimum height
    of the inlet patch; therefore, the minimum of \c z is always zero
    irrespective of the absolute z-coordinate of the computational patch.
    - The derived ABL expressions automatically satisfy the simplified transport
    equation for \c k. Yet the same expressions only satisfy the simplified
    transport equation for \c epsilon when the model constants \c sigmaEpsilon
    is 1.11 with \c kappa=0.4 (\c HW:p. 358).
    - \c atmBoundaryLayerInlet boundary conditions inherit \c inletOutlet
    traits, so that a given inlet condition can be supplied from all sides of
    the domain, e.g. a ground-normal cylinder domain having a single
    inlet/outlet boundary where the changes between inlet and outlet depend
    on the wind direction and patch normals, so that any change in inflow
    orientation can be handled with the same mesh.
    - \c d is the displacement height, and "is relevant for flows over forests
    and cities" (E:p. 28). "The displacement height gives the vertical
    displacement of the entire flow regime over areas which are densely covered
    with obstacles such as trees or buildings" (E:p. 28).

SourceFiles
    atmBoundaryLayerYGJ.C

\*---------------------------------------------------------------------------*/

#ifndef atmBoundaryLayerYGJ_H
#define atmBoundaryLayerYGJ_H

#include "fvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class atmBoundaryLayerYGJ Declaration
\*---------------------------------------------------------------------------*/

class atmBoundaryLayerYGJ
{
    // Private static data

        //- Default value of the Von Karman constant
        static const scalar kappaDefault_;

        //- Default value of the turbulent viscosity coefficient
        static const scalar CmuDefault_;

        //- Default value of the curve-fitting coefficient
        static const scalar C1Default_;

        //- Default value of the curve-fitting coefficient
        static const scalar C2Default_;

    // Private Data

        //- Flow direction
        vector flowDir_;

        //- Direction of the z-coordinate
        vector zDir_;

        //- Von Karman constant
        const scalar kappa_;

        //- Turbulent viscosity coefficient
        const scalar Cmu_;

        //- Curve-fitting coefficient
        const scalar C1_;

        //- Curve-fitting coefficient
        const scalar C2_;

        //- Reference velocity
        const scalar Uref_;

        //- Reference height
        const scalar Zref_;

        //- Surface roughness height
        scalarField z0_;

        //- Minimum coordinate value in z direction
        scalarField zGround_;

        //- Friction velocity
        scalarField Ustar_;

        //- True if the boundary layer is offset within the domain
        bool offset_;

        //- Velocity of the flow below the boundary layer
        const scalar Ulower_;

        //- Turbulence kinetic energy of the flow below the boundary layer
        const scalar kLower_;

        //- Turbulence kinetic energy dissipation rate of the flow below the
        //  boundary layer
        const scalar epsilonLower_;


    // Private Member Functions

        //- Initialisation shared by multiple constructors
        void init();


public:

    // Constructors

        //- Construct null
        atmBoundaryLayerYGJ();

        //- Construct from components
        atmBoundaryLayerYGJ
        (
            const vector& flowDir,
            const vector& zDir,
            const scalar Uref,
            const scalar Zref,
            const scalarField& z0,
            const scalarField& zGround,
            const scalar kappa = kappaDefault_,
            const scalar Cmu = CmuDefault_,
            const scalar C1 = C1Default_,
            const scalar C2 = C2Default_,
            const scalar ULower = 0,
            const scalar kLower = 0,
            const scalar epsilonLower = 0
        );

        //- Construct from the coordinates field and dictionary
        atmBoundaryLayerYGJ(const vectorField& p, const dictionary&);

        //- Construct by mapping given
        // atmBoundaryLayerYGJ onto a new patch
        atmBoundaryLayerYGJ
        (
            const atmBoundaryLayerYGJ&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmBoundaryLayerYGJ(const atmBoundaryLayerYGJ&);


    // Member Functions

        // Access

            //- Return flow direction
            const vector& flowDir() const
            {
                return flowDir_;
            }

            //- Return the ground-normal direction
            const vector& zDir() const
            {
                return zDir_;
            }

            //- Return friction velocity
            const scalarField& Ustar() const
            {
                return Ustar_;
            }


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            void rmap(const atmBoundaryLayerYGJ&, const labelList&);


        // Evaluation

            //- Return the velocity distribution for the ATM
            tmp<vectorField> U(const vectorField& p) const;

            //- Return the turbulent kinetic energy distribution for the ATM
            tmp<scalarField> k(const vectorField& p) const;

            //- Return the turbulent dissipation rate distribution for the ATM
            tmp<scalarField> epsilon(const vectorField& p) const;


        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
