/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::AbInitio

Description
     Ab initio (AI) potentials baed BinaryCollision Model

\*---------------------------------------------------------------------------*/

#ifndef AbInitio_H
#define AbInitio_H

#include "BinaryCollisionModel.H"
#include "dsmcCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class AbInitio Declaration
\*---------------------------------------------------------------------------*/


class AbInitio
:
    public BinaryCollisionModel
{
    // Private data
    //
        class AbInitioMatrix :
        {
            // priviate
            const word deflectionAngleTableFileName_;
            const label numRows_;
            const label numColumns_;
            const scalar G_;
            const char separator_;
        public:
            AbInitioMatrix(const IOdictionary& dict);
        }

        //- The coefficients dictionary
        const dictionary coeffDict_;

        List<word> deflectionAngleTableFileNames_;

        // store the deflection angle table 2D
        scalarField deflectionAngleTable_;
        // store the total cross section angle 1D
        scalarField sigmaTtable_;
        label numRows_;
        label numColumns_;
        scalar G_;
        const char separator_;

public:

    //- Runtime type information
    TypeName("AbInitio");


    // Constructors

        //- Construct from dictionary
        AbInitio
        (
            const dictionary& dict,
            dsmcCloud& cloud
        );


    //- Destructor
    ~AbInitio();


    // Member Functions

        //- Flag to indicate whether model activates collision model
        bool active() const;

        //- Return the collision cross section * relative velocity product
        scalar sigmaTcR
        (
            const dsmcParcel& pP,
            const dsmcParcel& pQ
        ) const;

        //- Apply collision
        void collide
        (
            dsmcParcel& pP,
            dsmcParcel& pQ,
            label& cellI
        );

        const dictionary& coeffDict() const;

        label deflectionAngleRow(scalar cR) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "AbInitio.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
