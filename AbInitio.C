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

\*---------------------------------------------------------------------------*/

#include "AbInitio.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(AbInitio, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, AbInitio, dictionary);
};


Foam::AbInitio::AbInitioMatrix
(
    const dictionary& dict
)
:


Foam::AbInitio::AbInitio
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_(dict.subDict(typeName + "Coeffs")),
    deflectionAngleTableFileName_(coeffDict_.lookup("deflectionAngleTableFileName")),
    numRows_(readLabel(coeffDict_.lookup("numRows"))),
    numColumns_(readLabel(coeffDict_.lookup("numColumns"))),
    G_(readScalar(coeffDict_.lookup("G"))),
    separator_(coeffDict_.lookupOrDefault<string>("separator",string(","))[0])
{
    // allocate the memory
    deflectionAngleTable_.setSize(numRows_*numColumns_);
    sigmaTtable_.setSize(numRows_);

    // read the tables
    IFstream in(cloud.time().constant()/deflectionAngleTableFileName_);

    // records the current line in the table
    label lines = 0;
    while (in.good())
    {
        string line;
        in.getLine(line);

        DynamicList<string> splitted;

        std::size_t pos = 0;
        while (pos != std::string::npos)
        {
            std::size_t nPos = line.find(separator_, pos);

            if (nPos == std::string::npos)
            {
                splitted.append(line.substr(pos));
                pos=nPos;
            }
            else
            {
                splitted.append(line.substr(pos, nPos-pos));
                pos=nPos+1;
            }
        }

        if (splitted.size() <= 1)
        {
            break;
        }

        // first numColumns are deflection angles
        for (label i=0; i < numColumns_; i++)
        {
            deflectionAngleTable_[lines*numColumns_ + i] =
                 readScalar(IStringStream(splitted[i])());
        }

        // last by 1 column is deflection angles
        sigmaTtable_[lines] = readScalar(IStringStream(splitted[numColumns_])())*1e-20;

        // next line ready
        lines += 1;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::AbInitio::~AbInitio()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::AbInitio::active() const
{
    return true;
}



Foam::scalar Foam::AbInitio::sigmaTcR
(
    const dsmcParcel& pP,
    const dsmcParcel& pQ
) const
{
//     const dmscCloud& cloud(this->owner());


//     if(typeIdP == 1 && typeIdQ == 0)
//     {
//         omegaPQ = 0.725;
//     }
//     
//     if(typeIdP == 0 && typeIdQ == 1)
//     {
//         omegaPQ = 0.725;
//     }


    scalar cR = mag(pP.U() - pQ.U());

    if (cR < VSMALL)
    {
        return 0;
    }

    // lookup table
    label row = deflectionAngleRow(cR); 
    scalar sigmaTPQ = sigmaTtable_[row];
    return sigmaTPQ*cR;
}

void Foam::AbInitio::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    label& cellI
)
{
//     dsmcCloud& cloud_(this->owner());

    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    
//     if(typeIdP == 1 && typeIdQ == 0)
//     {
//         alphaPQ = 1.0/1.64;
//     }
//     
//     if(typeIdP == 0 && typeIdQ == 1)
//     {
//         alphaPQ = 1.0/1.64;
//     }
    
//     Info << "alphaPQ = " << alphaPQ << endl;
        
    scalar collisionSeparation = sqrt(
            sqr(pP.position().x() - pQ.position().x()) +
            sqr(pP.position().y() - pQ.position().y())
    );
    
    cloud_.cellPropMeasurements().collisionSeparation()[cellI] += 
                                                        collisionSeparation;
    cloud_.cellPropMeasurements().nColls()[cellI]++;

    Random& rndGen(cloud_.rndGen());

    scalar mP = cloud_.constProps(typeIdP).mass();

    scalar mQ = cloud_.constProps(typeIdQ).mass();

    vector Ucm = (mP*UP + mQ*UQ)/(mP + mQ);

    scalar cR = mag(UP - UQ);
    
    vector cRComponents = UP - UQ;

    label row = deflectionAngleRow(cR); 
    label col = rndGen.integer(1,numColumns_);
    // lookup table
    scalar cosTheta = deflectionAngleTable_[row*numColumns_ + col - 1];

// lhzhu: we needs to lookup theta from the pre-stored matrix based on cR
    // store the matrix as a run time constructed properties of the collision model

    //scalar cosTheta = 2.0*(pow(rndGen.scalar01(),(1.0/alphaPQ))) - 1.0;
    scalar sinTheta = sqrt(1.0 -cosTheta*cosTheta);

    scalar phi = twoPi*rndGen.scalar01();
    
    scalar D = sqrt(cRComponents.y()*cRComponents.y() + cRComponents.z()*cRComponents.z());
    
    vector postCollisionRelU = vector::zero;
    
//     if(D > VSMALL)
//     {
        postCollisionRelU =
            vector
            (
                cosTheta*cRComponents.x() + sinTheta*sin(phi)*D,
                cosTheta*cRComponents.y() + sinTheta*(cR*cRComponents.z()*cos(phi) - cRComponents.x()*cRComponents.y()*sin(phi))/D,
                cosTheta*cRComponents.z() - sinTheta*(cR*cRComponents.y()*cos(phi) + cRComponents.x()*cRComponents.z()*sin(phi))/D
            ); //Bird, equation 2.22
//     }
//     else
//     {
//         postCollisionRelU =
//             vector
//             (
//                 cosTheta*cRComponents.x(),
//                 sinTheta*cos(phi)*cRComponents.x(),
//                 sinTheta*sin(phi)*cRComponents.x()
//             );
//     }

    UP = Ucm + postCollisionRelU*mQ/(mP + mQ);

    UQ = Ucm - postCollisionRelU*mP/(mP + mQ);
    
    label classificationP = pP.classification();
    label classificationQ = pQ.classification();
    
    //- class I molecule changes to class
    //- III molecule when it collides with either class II or class III
    //- molecules.
    
    if(classificationP == 0 && classificationQ == 1)
    {
        pP.classification() = 2;
    }
    
    if(classificationQ == 0 && classificationP == 1)
    {
        pQ.classification() = 2;
    }
    
    if(classificationP == 0 && classificationQ == 2)
    {
        pP.classification() = 2;
    }
    
    if(classificationQ == 0 && classificationP == 2)
    {
        pQ.classification() = 2;
    }
}

const Foam::dictionary& Foam::AbInitio::coeffDict() const
{
    return coeffDict_;
}
// ************************************************************************* //
//

Foam::label Foam::AbInitio::deflectionAngleRow(Foam::scalar cR) const
{
    Foam::label j = floor(log(1+cR/G_)/log(1.005) + 0.5);
    if (j>numRows_) j = numRows_;
    return j - 1;
}
