/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "oka.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::oka<CloudType>::applyToPatch
(
    const label globalPatchi
) const
{
    forAll(patchIDs_, i)
    {
        if (patchIDs_[i] == globalPatchi)
        {
            return i;
        }
    }

    return -1;
}


template<class CloudType>
void Foam::oka<CloudType>::write()
{
    if (QPtr_.valid())
    {
        QPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "QPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::oka<CloudType>::oka
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    QPtr_(nullptr),
    patchIDs_(),
    p_(readScalar(this->coeffDict().lookup("p"))),
    psi_(this->coeffDict().template lookupOrDefault<scalar>("psi", 2.0)),
    K_(this->coeffDict().template lookupOrDefault<scalar>("K", 2.0))
{
    const wordList allPatchNames(owner.mesh().boundaryMesh().names());
    const wordReList patchNames(this->coeffDict().lookup("patches"));

    labelHashSet uniqIds;
    for (const wordRe& re : patchNames)
    {
        labelList ids = findStrings(re, allPatchNames);

        if (ids.empty())
        {
            WarningInFunction
                << "Cannot find any patch names matching " << re
                << endl;
        }

        uniqIds.insert(ids);
    }

    patchIDs_ = uniqIds.sortedToc();

    // trigger creation of the Q field
    preEvolve();
}


template<class CloudType>
Foam::oka<CloudType>::oka
(
    const oka<CloudType>& pe
)
:
    CloudFunctionObject<CloudType>(pe),
    QPtr_(nullptr),
    patchIDs_(pe.patchIDs_),
    p_(pe.p_),
    psi_(pe.psi_),
    K_(pe.K_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::oka<CloudType>::preEvolve()
{
    if (QPtr_.valid())
    {
        QPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        QPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Q",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimVolume, Zero)
            )
        );
    }
}


template<class CloudType>
void Foam::oka<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    bool&
)
{
    const label patchi = pp.index();

    const label localPatchi = applyToPatch(patchi);

    if (localPatchi != -1)
    {
        vector nw;
        vector Up;

        // patch-normal direction
        this->owner().patchData(p, pp, nw, Up);

        // particle velocity relative to patch
        const vector& U = p.U() - Up;

        // quick reject if particle travelling away from the patch
        if ((nw & U) < 0)
        {
            return;
        }

        const scalar magU = mag(U);
        const vector Udir = U/magU;

        // determine impact angle, alpha
        const scalar alpha = mathematical::pi/2.0 - acos(nw & Udir);
	
         const scalar coeff =p.nParticle()*(pow(sin(alpha),0.62))*(pow(1+(0.4*(1-sin(alpha))),5.68));//*sqr(magU)/(p_*psi_*K_);//const scalar coeff = 
       // const scalar coeff = p.nParticle()*p.mass()*sqr(magU)/(p_*psi_*K_);
	
        const label patchFacei = pp.whichFace(p.face());
	//const scalarField magSf(mag(patchFace.faceAreas()));
	//scalar ar =       pp.magSf();//[patchFacei];
        scalar& Q = QPtr_->boundaryFieldRef()[patchi][patchFacei];
        //scalar ern =65*(pow((0.019*0.4),(-0.12*1.58)))*pow((magU/104),2.2)*pow((10/326),0.19);
        const fvMesh& mesh = this->owner().mesh();
	//const polyPatch& patch = mesh.boundaryMesh()[patchi];
    //const pointField& points = patch.points();
	//const scalarField magSf(mag(patch.faceAreas()));
	const surfaceScalarField& magSf = mesh.magSf();//scalar arsu =magSf[patchFacei];//const surfaceScalarField& magSf = mesh.magSf(); scalar area = sum(magSf[patchi]);
	scalar area=(magSf[patchFacei]);//[patchFacei]);
            Q += ((coeff*(2.814399483e-12)*pow(magU,2.2))/area);///(area*1.413716694e-7));//65*(pow((0.019*0.4),(-0.12*1.58)))*pow((magU/104),2.2)*pow((10/326),0.19)*(pow(sin(alpha),0.62))*(pow(1+(0.4*(1-sin(alpha))),5.68));//coeff*ern*gal;
        
        
    }
}


// ************************************************************************* //
