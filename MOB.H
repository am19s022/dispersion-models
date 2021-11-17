/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::MOB

Description
    The velocity is perturbed in the direction of -grad(k), with a
    Gaussian random number distribution with variance sigma.
    where sigma is defined below

\*---------------------------------------------------------------------------*/

#ifndef MOB_H
#define MOB_H

#include "DispersionRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class MOB Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MOB
:
    public DispersionRASModel<CloudType>
{
protected:

    // Protected data

        // Locally cached turbulence fields
		volScalarField SumPar_;
		volVectorField UPar_;
            //- Gradient of k
	    const volVectorField* gradkPtr_;
            const volTensorField* rsPtr_;
	    const volVectorField* divRsPtr_;
	    //const volVectorField* gradkPtr_;
            //- Take ownership of the grad(k)
            mutable bool ownrs_;
	     mutable bool owndivRs_;
	mutable bool ownGradK_;
          const volScalarField* diaPtr_;
	 mutable bool owndia_;  // mutable bool ownGradk_;
	const volVectorField* UensPtr_;
	mutable bool ownUens_;
	const volScalarField* spPtr_;
	 mutable bool ownsp_;

public:

    //- Runtime type information
    TypeName("MOB");


    // Constructors

        //- Construct from components
        MOB(const dictionary& dict, CloudType& owner);

        //- Construct copy
        MOB(const MOB<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const
        {
            return autoPtr<DispersionModel<CloudType>>
            (
                new MOB<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MOB();


    // Member Functions

        //- Cache carrier fields
        virtual void cacheFields(const bool store);

        //- Update (disperse particles)
        virtual vector update
        (
            const scalar dt,
            const label celli,
            const vector& U,
            const vector& Uc,
            vector& UTurb,
            scalar& tTurb
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MOB.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
