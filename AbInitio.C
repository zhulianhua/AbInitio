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
#include "Switch.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


namespace Foam
{
    defineTypeNameAndDebug(AbInitio, 0);
    addToRunTimeSelectionTable(BinaryCollisionModel, AbInitio, dictionary);
};


Foam::AbInitio::AbInitioMatrix::AbInitioMatrix
(
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    cloud_(cloud),
    deflectionAngleCosinTableFileName_(dict.lookup("deflectionAngleCosinTableFileName")),
    numRows_(readLabel(dict.lookup("numRows"))),
    numColumns_(readLabel(dict.lookup("numColumns"))),
    G_(readScalar(dict.lookup("G"))),
    separator_(dict.lookupOrDefault<string>("separator",string(","))[0]),
    useInterpolatedSigmaT_(dict.lookupOrDefault<Switch>("useInterpolatedSigmaT", false))
{
    // allocate the memory
    deflectionAngleCosinTable_.setSize(numRows_*numColumns_);
    sigmaTtable_.setSize(numRows_);

    // read the tables
    IFstream in(cloud.time().constant()/deflectionAngleCosinTableFileName_);

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
            deflectionAngleCosinTable_[lines*numColumns_ + i] =
                 readScalar(IStringStream(splitted[i])());
        }

        // last by 1 column is total cross section in unit of 1e-20 m^2
        sigmaTtable_[lines] = readScalar(IStringStream(splitted[numColumns_])())*1e-20;

        // next line ready
        lines += 1;
    }
}

inline Foam::scalar Foam::AbInitio::AbInitioMatrix::deflectionAngleCosin ( Foam::scalar cR) const
{
    Random& rndGen(cloud_.rndGen());
    label row = deflectionAngleRow(cR);
    label col = rndGen.integer(1,numColumns_);
    // lookup table
    return deflectionAngleCosinTable_[row*numColumns_ + col - 1];
}

inline Foam::scalar Foam::AbInitio::AbInitioMatrix::sigmaT ( Foam::scalar cR) const
{
    if(useInterpolatedSigmaT_) {
        label rowA, rowB;
        deflectionAngleRowBetween(cR, rowA, rowB);
        // recover the original cRA and cRB
        // NOTE: can be pre-stored also in that matrix
        scalar cRA = G_*(pow(1.005, rowA+1)-1.0);
        scalar cRB = G_*(pow(1.005, rowB+1)-1.0);
        return (cR - cRB)*(sigmaTtable_[rowA] - sigmaTtable_[rowB])/(cRA - cRB)
            + sigmaTtable_[rowB];
    }
    else{
        label row = deflectionAngleRow(cR);
        return sigmaTtable_[row];
    }
}

inline Foam::label Foam::AbInitio::AbInitioMatrix::deflectionAngleRow(Foam::scalar cR) const
{
    Foam::label j = floor(log(1+cR/G_)/log(1.005) + 0.5);
    if (j>numRows_) j = numRows_;
    return j - 1;
}

void Foam::AbInitio::AbInitioMatrix::deflectionAngleRowBetween(Foam::scalar cR, Foam::label& rowA, Foam::label& rowB) const
{
    Foam::label jb = floor(log(1+cR/G_)/log(1.005) + 0.5);
    Foam::label ja = jb+1;

    if (jb>numRows_) 
    {
        ja = numRows_;
        jb = numRows_;
    }

    rowA =  ja-1;
    rowB =  jb-1;
}

Foam::AbInitio::AbInitio
(
    const dictionary& dict,
    dsmcCloud& cloud
)
:
    BinaryCollisionModel(dict, cloud),
    coeffDict_(dict.subDict(typeName + "Coeffs"))
{
    label nComponents = cloud.typeIdList().size();
    // set the list size
    AImatrixs_.setSize(nComponents);

    for(label i = 0; i < nComponents; i++)
    {
        AImatrixs_.set(i, new PtrList<AbInitioMatrix>);
    }

    // set the sub list size
    for(label i = 0; i < nComponents; i++)
    {
        AImatrixs_[i].setSize(nComponents);
    }

    // construct the matrixs
    for(label j = 0; j < nComponents; j++)
    {
        for(label i = j; i < nComponents; i++)
        {
            word componentP = cloud.typeIdList()[j];
            word componentQ = cloud.typeIdList()[i];
            AImatrixs_[j].set(i, new AbInitioMatrix(cloud, coeffDict_.subDict(componentP+'-'+componentQ)));
            // the symmetry list element
            if (i > j) AImatrixs_[i].set(j, new AbInitioMatrix(cloud, coeffDict_.subDict(componentP+'-'+componentQ)));
        }
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
    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();

    scalar cR = mag(pP.U() - pQ.U());

    return AImatrixs_[typeIdP][typeIdQ].sigmaT(cR)*cR;
}

void Foam::AbInitio::collide
(
    dsmcParcel& pP,
    dsmcParcel& pQ,
    label& cellI
)
{
    label typeIdP = pP.typeId();
    label typeIdQ = pQ.typeId();
    vector& UP = pP.U();
    vector& UQ = pQ.U();
    
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

    scalar cosTheta = AImatrixs_[typeIdP][typeIdQ].deflectionAngleCosin(cR);
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

const Foam::AbInitio::AbInitioMatrix& Foam::AbInitio::AImatrixs(label j, label i) const
{
    return AImatrixs_[j][i];
}
