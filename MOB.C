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

\*---------------------------------------------------------------------------*/

#include "MOB.H"
#include "demandDrivenData.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "KinematicParcel.H"
#include "ParticleErosion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// this is CRW in name of MOB for making
template<class CloudType>
Foam::MOB<CloudType>::MOB
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner),
    SumPar_
    (
        IOobject
        (
            this->owner().name() + ":SumPar",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->owner().mesh(),
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),	
	UPar_
    (
        IOobject
        (
            this->owner().name() + ":UPar",
            this->owner().db().time().timeName(),
            this->owner().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->owner().mesh(),
        dimensionedVector(dimVelocity, Zero),
        zeroGradientFvPatchVectorField::typeName
    ),
    rsPtr_(nullptr),
    ownrs_(false),
    divRsPtr_(nullptr),
    owndivRs_(false),
    gradkPtr_(nullptr),
    ownGradK_(false),
	diaPtr_(nullptr),
    owndia_(false),
    UensPtr_(nullptr),
    ownUens_(false),
	spPtr_(nullptr),
    ownsp_(false)                           //gradkPtr_(nullptr),
                               // ownGradK_(false)
{SumPar_ = this->owner().SumPar();
    SumPar_.oldTime();
	UPar_ = this->owner().UPar();
    UPar_.oldTime();
}


template<class CloudType>
Foam::MOB<CloudType>::MOB
(
    const MOB<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm),
	SumPar_(dm.SumPar_),
	UPar_(dm.UPar_),
    rsPtr_(dm.rsPtr_),
    ownrs_(dm.ownrs_),
    divRsPtr_(dm.divRsPtr_),
    owndivRs_(dm.owndivRs_),
    gradkPtr_(dm.gradkPtr_),
    ownGradK_(dm.ownGradK_),
   // gradkPtr_(dm.gradkPtr_),
    //ownGradK_(dm.ownG#include "fvmDiv.H"radK_)
    diaPtr_(dm.diaPtr_),
    owndia_(dm.owndia_),
    UensPtr_(dm.UensPtr_),
    ownUens_(dm.ownUens_),
	spPtr_(dm.spPtr_),
    ownsp_(dm.ownsp_)

{     SumPar_.oldTime();
	UPar_.oldTime();
   dm.owndivRs_ = false;
    dm.ownGradK_ = false;// dm.ownGradK_ = false;
    dm.ownrs_ = false;
	dm.owndia_ = false;
   dm.ownUens_ = false;
	dm.ownsp_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MOB<CloudType>::~MOB()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MOB<CloudType>::cacheFields(const bool store)
{
    DispersionRASModel<CloudType>::cacheFields(store);

    if (store)
    {
       rsPtr_ = ((*this->nutPtr_)*(  -(fvc::grad(*this->UPtr_))  -  (T(fvc::grad(*this->UPtr_)) )) ).ptr();//+((2/3)*I*(*this->kPtr_))).ptr(); // MAIN STEP
        ownrs_ = true;
       // RsPtr_=rsPtr_;
	divRsPtr_ = fvc::div(*this->rsPtr_).ptr(); // MAIN STEP
        owndivRs_ = true;
        gradkPtr_ = fvc::grad(*this->kPtr_).ptr();
        ownGradK_ = true;
        //gradkPtr_ = fvc::grad(*this->kPtr_).ptr();
        //ownGradK_ = true;
	diaPtr_ = (this->owner().dia()).ptr();
	owndia_ = true;
	UensPtr_ = (UPar_ + (this->owner().UPar())).ptr();
        ownUens_ = true;
	spPtr_ = (SumPar_ + (this->owner().SumPar())).ptr();
        ownsp_ = true;
    }
    else
    {
        if (ownrs_)
        {
            deleteDemandDrivenData(rsPtr_);
            rsPtr_ = nullptr;
            ownrs_ = false;
	   
        }

	if (owndivRs_)
        {
            deleteDemandDrivenData(divRsPtr_);
            divRsPtr_ = nullptr;
            owndivRs_ = false;
	   
        }
        
        if (ownGradK_)
        {
            deleteDemandDrivenData(gradkPtr_);
            gradkPtr_ = nullptr;
            ownGradK_ = false;
        }
          if (owndia_)
         {
            deleteDemandDrivenData(diaPtr_);
            diaPtr_ = nullptr;
            owndia_ = false;
        }
	if (ownUens_)
        {
            deleteDemandDrivenData(UensPtr_);
            UensPtr_ = nullptr;
            ownUens_ = false;
        }
	if (ownsp_)
        {
            deleteDemandDrivenData(spPtr_);
            spPtr_ = nullptr;
            ownsp_ = false;
        }

    }

}


template<class CloudType>
Foam::vector Foam::MOB<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    scalar& tTurb
)
{
    Random& rnd = this->owner().rndGen();
    //volScalarField ptot=(this->owner().SumPar())();
     //volScalarField ptot1 =ptot.oldTime() ;
	//volScalarField ptot2=ptot + ptot1;
	//const fvMesh& mesh = this->owner().mesh();
	//tensor mn =fvc::grad(Uc) ;       							
    const scalar cps = 0.16432;
    const vector& Uens = this->UensPtr_->primitiveField()[celli];
	const scalar& sp = this->spPtr_->primitiveField()[celli];
    vector Utem= Uens;scalar sptem=sp;
	if(sptem>0)
       {
            Utem = Utem/(sptem+ROOTVSMALL);
        } 
    const scalar k = this->kPtr_->primitiveField()[celli];
	const scalar nut = this->nutPtr_->primitiveField()[celli];
    const scalar epsilon =
        this->epsilonPtr_->primitiveField()[celli] + ROOTVSMALL;
	
    tensor rs = this->rsPtr_->primitiveField()[celli];
    tensor ts = rs;//check if nut or mut
    ts.xx()=ts.xx()+(2.0*k/3.0);//restr calculation
    ts.yy()=ts.yy()+(2.0*k/3.0);
    ts.zz()=ts.zz()+(2.0*k/3.0);
    const vector& divRs = this->divRsPtr_->primitiveField()[celli];
    const vector& gradk = this->gradkPtr_->primitiveField()[celli];
    vector drs = divRs;
    drs = drs +((2/3)*gradk) ;
    const scalar UrelMag = mag(U - Uc - UTurb);
    const scalar tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + SMALL));//UTurb = Uc/100;
    vector Ur = -Utem + Uc + UTurb ;//formation of new coordinate system
    vector Ur1=-Ur;
    scalar ctet = Ur1.x()/(ROOTVSMALL+sqrt(sqr(Ur1.x())+sqr(Ur1.y())));
    scalar stet= -Ur1.y()/(ROOTVSMALL+sqrt(sqr(Ur1.x())+sqr(Ur1.y())));
    scalar cph= sqrt(sqr(Ur1.x())+sqr(Ur1.y()))/(ROOTVSMALL+mag(Ur1));	
    scalar sph = - Ur1.z()/(ROOTVSMALL+mag(Ur1));
    tensor C(ctet*cph,stet,ctet*sph,-stet*cph,ctet,-stet*sph,-sph,0,cph);C = C.T();
    tensor rtr = transform(C,ts) ;
    scalar taul = 0.4819*k/epsilon ;//taul=tTurbLoc;
    scalar tp =2700*(sqr(10e-6))/(18*1.2*1.57e-5);
    scalar rerel = (1e-5)*UrelMag/(1.57e-5);
    scalar coe=0;
	if (rerel<=1)
	{coe =1 ;}
	else if (rerel<=1000)
	{coe=1+(0.15*pow(rerel,0.687));}
	else
	{coe=0.44*rerel/24 ;}
    tp =tp/coe;
    scalar ste = 0.356*tp/(ROOTVSMALL+taul) ;
    scalar ste1=pow((1+ste),(-0.4-(0.004*ste)));
    scalar bs=1-(0.644*ste1);
    scalar erisq=sqr(mag(Ur))/((2*k/3)+ROOTVSMALL);
    scalar eri = sqrt(erisq);
    scalar den1 = sqrt(1+sqr(bs*eri));
    scalar den2 = sqrt(1+sqr(2*bs*eri)); 
    scalar taus1 = (bs/0.356)*taul/(ROOTVSMALL+den1); 
    scalar taus2 = (bs/0.356)*taul/(ROOTVSMALL+den2);
    //dt to be hard coded as 5e-6  
	scalar  a1=exp(-5e-6/(taus1+ROOTVSMALL));//(ROOTVSMALL+taus1));
   scalar a2=exp(-5e-6/(taus2+ROOTVSMALL));//ROOTVSMALL+taus2));
    scalar a3 =a2;
    tensor aij(a1,0,0,0,a2,0,0,0,a3);	
    tensor diff(0,0,0,0,0,0,0,0,0);	
    /*diff.xx()=(1-(aij.xx()*aij.xx()))*rtr.xx();
    diff.xy()=(1-(aij.xx()*aij.yy()))*rtr.xy();
    diff.xz()=(1-(aij.xx()*aij.zz()))*rtr.xz();              
    diff.yx()=(1-(aij.yy()*aij.xx()))*rtr.yx();
    diff.yy()=(1-(aij.yy()*aij.yy()))*rtr.yy();
    diff.yz()=(1-(aij.yy()*aij.zz()))*rtr.yz();
    diff.zx()=(1-(aij.zz()*aij.xx()))*rtr.zx();
    diff.zy()=(1-(aij.zz()*aij.yy()))*rtr.zy();
    diff.zz()=(1-(aij.zz()*aij.zz()))*rtr.zz();*/
    
    //tensor  diff=rtr;
    diff.xx()=(1-(a1*a1))*rtr.xx();
    diff.xy()=(1-(a1*a2))*rtr.xy();
    diff.xz()=(1-(a1*a3))*rtr.xz();	
    diff.yx()=(1-(a2*a1))*rtr.yx();   
    diff.yy()=(1-(a2*a2))*rtr.yy(); 
    diff.yz()=(1-(a2*a3))*rtr.yz();
    diff.zx()=(1-(a3*a1))*rtr.zx();
    diff.zy()=(1-(a3*a2))*rtr.zy();
    diff.zz()=(1-(a3*a3))*rtr.zz();
    tensor l(0,0,0,0,0,0,0,0,0);
    l.xx()=sqrt(diff.xx());
    l.xy()=0;     
    l.xz()=0;
    l.yx()=diff.yx()/(ROOTVSMALL+l.xx());
    l.yy()=sqrt(diff.yy()-sqr(l.yx()));
    l.yz()=0;
    l.zx()=diff.zx()/l.xx() ;
    l.zy() = ( diff.zy() - (l.zx()*l.yx()))/(l.yy()+ROOTVSMALL);
    l.zz()= sqrt(diff.zz()-sqr(l.zx())-sqr(l.zy()));	
 
//crw end	
	

     scalar sno=tp/((k/(7*epsilon))+ROOTVSMALL);//drift correction
	vector dcr = drs/(1+sno);
	vector Utr = transform(C,UTurb);       
	scalar cz1 = rnd.GaussNormal<scalar>(); scalar cz2 = rnd.GaussNormal<scalar>(); scalar cz3 = rnd.GaussNormal<scalar>();
	vector zi(0,0,0); zi.x() = cz1; zi.y() = cz2; zi.z() = cz3;
	vector bi(0,0,0); bi.x()= (l.xx()*cz1);bi.y()=(l.yx()*cz1)+(l.yy()*cz2);bi.z()=(l.zx()*cz1)+(l.zy()*cz2)+(l.zz()*cz3); 
	UTurb = (Utr & aij ) + bi;//(zi & l);
        UTurb = transform(C.T(),UTurb);
	UTurb = UTurb +(dt*dcr);
	return Uc + UTurb;



	
	
	/*vector Urp = U + Uc + UTurb ; 
 	vector x1(1,1,1); x1.x() = Ur.x()/(mag(Ur)+ROOTVSMALL); x1.y() = Ur.y()/(mag(Ur)+ROOTVSMALL); x1.y() = Ur.y()/(mag(Ur)+ROOTVSMALL);
	vector x2(1,1,1); x2.x() = Urp.x()/(mag(Urp)+ROOTVSMALL); x2.y() = Urp.y()/(mag(Urp)+ROOTVSMALL); x2.z()= Urp.z()/(mag(Urp)+ROOTVSMALL);
	vector per = Ur ^ Urp;
	vector x3(1,1,1); x3.x() = per.x()/(mag(per)+ROOTVSMALL); x3.y() = per.y()/(mag(per)+ROOTVSMALL); x3.z() = per.z()/(mag(per)+ROOTVSMALL);
	//error passage*/
	
	//isotropic crw

	//tensor C(x1.x(),x1.y(),x1.z(),x2.x(),x2.y(),x2.z(),x3.x(),x3.y(),x3.z()); 
	/*vector tau(0,0,0); tau.x() = 0.4819*k/epsilon; tau.y() = 0.4819*k/epsilon; tau.y() = 0.4819*k/epsilon;
	vector el(0,0,0); el.x() = sqrt(2*0.4819/(3*0.356))*pow(k,1.5)/epsilon; el.y() = 0.5*sqrt(2*0.4819/(3*0.356))*pow(k,1.5)/epsilon;
	el.z() = 0.5*sqrt(2*0.4819/(3*0.356))*pow(k,1.5)/epsilon;
	tensor aij(0,0,0,0,0,0,0,0,0); 
	aij.xx() = 
(-dt/tau.x())*exp(-dt*mag(Ur)/el.x())*(1-( (dt*mag(Ur)/(2*el.x()) ) * (1-sqr(Ur.x()/mag(Ur)))));
	aij.yy() = exp(-dt/tau.y())*exp(-dt*mag(Ur)/el.y())*(1-((dt*mag(Ur)/(2*el.y()) )*(1-sqr(Ur.y()/mag(Ur)))));
	aij.zz() = exp(-dt/tau.z())*exp(-dt*mag(Ur)/el.z())*(1-((dt*mag(Ur)/(2*el.z()) )*(1-sqr(Ur.z()/mag(Ur)))));
	tensor bij(0,0,0,0,0,0,0,0,0);	
	//tensor rtr = (C & ts ) & (C.T()) ;
	bij.xx() = sqrt(rtr.xx()*(1-(sqr(aij.xx()))));
	bij.yy() = sqrt(rtr.yy()*(1-(sqr(aij.yy()))));
	bij.zz() = sqrt(rtr.xx()*(1-(sqr(aij.xx()))));
	vector Utr = (C & UTurb);       
	scalar c1 = rnd.GaussNormal<scalar>(); scalar c2 = rnd.GaussNormal<scalar>(); scalar c3 = rnd.GaussNormal<scalar>();
	vector zi(0,0,0); zi.x() = c1; zi.y() = c2; zi.z() = c3;
	scalar tp = 2700*(sqr(10e-6))/(18*1.57e-5); scalar rerel = (2700e-5)*mag(Ur)/(1.57e-5); scalar coe=0;
	if (rerel<=1000)
	{coe=1+(0.15*pow(rerel,0.687));}
	else
	{coe=0.44*rerel/24 ;}
	tp =tp/coe;

        scalar sno = tp/((k/(7*epsilon))+ROOTVSMALL);//drift correction
	vector dcr = drs/(1+sno);
	//vector Utr = (C & UTurb);
	UTurb = (Utr & aij ) + (zi & bij);
	UTurb = ((C.T()) & Utr)+(dt*dcr);*/






	/*tensor rtr = (C & ts ) & (C.T()) ;// ( (T(C) & ts) & C );//transformation of re stress
	
           // until this checked

	vector ri(0,0,0);
	vector tai(0,0,0);tai.x()=0.4819*1.5*rtr.xx()/epsilon;tai.y()=0.4819*1.5*rtr.yy()/epsilon;tai.z()=0.4819*1.5*rtr.zz()/epsilon;
	

	//crw start
ri.x() = exp(-dt/(tai.x()+ROOTVSMALL) );                 //Rli calc
	ri.y() =  exp(-dt/(tai.y()+ROOTVSMALL) );
	ri.z()= exp(-dt/(tai.z()+ROOTVSMALL) );
	vector li(1,2,3);
	li.x()=(1.5*0.4819/(0.356))*((pow(rtr.xx(),1.5))/epsilon) ;
	li.y()=(1.5*0.4819/(0.5*0.356))*((pow(rtr.yy(),1.5))/epsilon);
	li.z()=(1.5*0.4819/(0.5*0.356))*((pow(rtr.zz(),1.5))/epsilon) ;  //Lei calc
	vector rrel = -Ur*dt;
	vector fi(1,2,3);
	fi.x()=exp(-mag(rrel)/(li.x()+ROOTVSMALL));
	fi.y()=exp(-mag(rrel)/(li.y()+ROOTVSMALL));
	fi.z()=exp(-mag(rrel)/(li.z()+ROOTVSMALL));
	vector gi(1,2,3);
	gi.x() = (1-(mag(rrel)/((2*li.x())+ROOTVSMALL)))*fi.x();	
	gi.y() = (1-(mag(rrel)/((2*li.y())+ROOTVSMALL)))*fi.y();	
	gi.z() = (1-(mag(rrel)/((2*li.z())+ROOTVSMALL)))*fi.z();
	vector rei(1,2,3);
	rei.x()=((fi.x()-gi.x())*(sqr(rrel.x()/(mag(rrel)+ROOTVSMALL)))) + gi.x() ;
	rei.y()=((fi.y()-gi.y())*(sqr(rrel.y()/(mag(rrel)+ROOTVSMALL)))) + gi.y() ;
	rei.z()=((fi.z()-gi.z())*(sqr(rrel.z()/(mag(rrel)+ROOTVSMALL)))) + gi.z() ;
	vector rsi(1,2,3);
	rsi.x() = ri.x()*rei.x();
	rsi.y() = ri.y()*rei.y();
	rsi.z() = ri.z()*rei.z();
	tensor aij(rsi.x(),0,0,0,rsi.y(),0,0,0,rsi.z());
	//tensor aijtr=(T(C)&(aij)&T);//for final re check
	tensor bij(1,0,0,0,1,0,0,0,1);
	bij.xx()=sqrt(rtr.xx()*(1-(sqr(rsi.x())) ));
	bij.yy()=sqrt(rtr.yy()*(1-(sqr(rsi.y()) ) ));	
	bij.zz()=sqrt(rtr.zz()*(1-(sqr(rsi.z()))));
	
//crw end	
	scalar tp =2700*(sqr(10e-6))/(18*1.57e-5);scalar rerel = (2700e-5)*mag(Ur)/(1.57e-5); scalar coe=0;
	if (rerel<=1000)
	{coe=1+(0.15*pow(rerel,0.687));}
	else
	{coe=0.44*rerel/24 ;}
	tp =tp/coe;

     scalar sno=tp/((k/(7*epsilon))+ROOTVSMALL);//drift correction
	vector dcr = drs/(1+sno);
//velocity transformations 
	vector Utr(0,0,0); Utr= C & UTurb ;scalar c1 =rnd.GaussNormal<scalar>();scalar c2= rnd.GaussNormal<scalar>();scalar c3=rnd.GaussNormal<scalar>();
	vector zi(0,0,0);zi.x()=c1;zi.y()=c2;zi.z()=c3;
	Utr = (Utr & aij ) + (zi & bij);Utr = (C.T()) & Utr;
	UTurb = Utr + (dt*drs);//& T(C);//Utr = Utr +(dt*dcr);UTurb= Utr;
	//UTurb = U/100;




//MOB start	
	/*vector sti(0,0,0);scalar tp =2700*(sqr(10e-6))/(18*1.57e-5);
	sti.x() = 0.356*tp/tai.x();//St ei calculation
	sti.y() = 2*0.356*tp/tai.y();
	sti.z() = 2*0.356*tp/tai.z();
	vector bs(0,0,0);//beta * calculation
	bs.x()= (0.356/0.356)*(tau.x()/tau.x())*(1-((1-(0.356*tau.x()/tau.x()))*pow(1+sti.x(),-0.4*(1+(0.01*sti.x())))));
	bs.y()= (2*0.356/0.356)*(tau.x()/tau.y())*(1-((1-(0.356*tau.y()/tau.x()))*pow(1+sti.x(),-0.4*(1+(0.01*sti.x())))));
	bs.x()= (2*0.356/0.356)*(tau.x()/tau.z())*(1-((1-(0.356*tau.z()/tau.x()))*pow(1+sti.x(),-0.4*(1+(0.01*sti.x())))));
	vector eri(0,0,0);//eri calc
	eri.x() = sqr(mag(Ur))/rtr.xx();
	eri.y() = sqr(mag(Ur))/rtr.yy();
	eri.z() = sqr(mag(Ur))/rtr.zz();
	vector tausi(0,0,0);//tau si calculation
	tausi.x() =(bs.x()/0.356)*tai.x()/sqrt(1+sqr((bs.x()*eri.x())));//tausi calculation
	tausi.y() =(bs.x()/0.356)*tai.y()/sqrt(1+sqr((bs.y()*eri.y())));
	tausi.z() =(bs.z()/0.356)*tai.z()/sqrt(1+sqr((bs.z()*eri.z())));
	vector ai(0,0,0)	;//ai calculation
	ai.x()=exp(-dt/tausi.x());// ai calculation
	ai.y()=exp(-dt/tausi.y());		
	ai.z()=exp(-dt/tausi.z());
	//choleksky
`	tensor cho(0,0,0,0,0,0,0,0,0);
	cho.xx() =(1-(ai.x()*ai.x()))*rtr.xx();
	cho.xy() =(1-(ai.x()*ai.y()))*rtr.xy();
	cho.xz() =(1-(ai.x()*ai.z()))*rtr.xz();
	cho.yx() =(1-(ai.y()*ai.x()))*rtr.yx();
	cho.yy() =(1-(ai.y()*ai.y()))*rtr.yy();
	cho.yz() =(1-(ai.y()*ai.z()))*rtr.yz();
	cho.zx() =(1-(ai.z()*ai.x()))*rtr.zx();
	cho.zy() =(1-(ai.z()*ai.y()))*rtr.zy();
	cho.zz() =(1-(ai.z()*ai.z()))*rtr.zz();
	//decomposition starts
	tensor low(0,0,0,0,0,0,0,0,0);
	low.xx()=sqrt(cho.xx());
	low.yx()=cho.yx()/low.xx();
	low.yy()=sqrt(cho.yy()-sqr(low.yx()));
	low.zx()=cho.zx()/low.xx();
	low.zy()=(cho.zy()-(cho.zx()*cho.yx()))/low.yy();
	low.yy()=sqrt(cho.zz()-sqr(low.zx())-sqr(low.zy()));
	vector bi(0,0,0);
	scalar c1 =rnd.GaussNormal<scalar>();scalar c2= rnd.GaussNormal<scalar>();scalar c3=rnd.GaussNormal<scalar>();
	vector zi(0,0,0);zi.x()=c1;zi.y()=c2;zi.z()=c3;
	bi = low & zi ;
	
	//drift correction term
	
	//final transformations
	//vector UTurbco= C & UTurb ;//term one forward transformation
	//vector vf(0,0,0);vf.x()=ai.x()*UTurbco.x();vf.y()=ai.y()*UTurbco.y();vf.z()=ai.z()*UTurbco.z();
	//vf =vf + (bi*dt) ;
	//vector vtrans = T(c)&vf;
	//vtrans=vtrans+dcr;//consideringdrift correction in global coordinates
//	UTurb=Uturb+vtrans;
//tensor rs
  //vector ih = divRs;ih = ih -((2/3)*gradk) ;
//ts = ts+(ts');
//vector rg= divRs;
//scalar gk = gradk;

 // Parcel is perturbed by the turbulence
   if (dt < tTurbLoc)
    {
        tTurb += dt;

        if (tTurb > tTurbLoc)
        {
            tTurb = 0.0;
	//ts.xx()=ts.xx()-(2.0*k/3.0);
	//ts.yy()=ts.yy()-(2.0*k/3.0);
	//ts.zz()=ts.zz()-(2.0*k/3.0);
	
	//vector v = fvc::grad(ts);
	//UTurb.x()= (-UTurb.x()*dt);/////+( (sqrt(2*ts.xx()))*rnd.GaussNormal<scalar>()) ;                               //sigma*fac*dir;
       	  // UTurb.y()= (-UTurb.y()*dt/1);//+( (sqrt(2*ts.yy()))*rnd.GaussNormal<scalar>()) ;
		// UTurb.z()= (-UTurb.z()*dt/1);///////+( (sqrt(2*ts.zz()/1))*rnd.GaussNormal<scalar>()) ;

            //const scalar sigma = sqrt(2.0*k/3.0);
            //const vector dir = -gradk/(mag(gradk) + SMALL);

           // scalar fac = 0.0;

            // In 2D calculations the -grad(k) is always
            // away from the axis of symmetry
            // This creates a 'hole' in the spray and to
            // prevent this we let fac be both negative/positive
            if (this->owner().mesh().nSolutionD() == 2)
            {
                fac = rnd.GaussNormal<scalar>();
            }
            else
            {
                fac = mag(rnd.GaussNormal<scalar>());
            }
	UTurb = U/100;
           
       }
    }
    else
    {
        tTurb = GREAT;
        UTurb = Zero;
    }*/
    //return Uc + UTurb;
}


// ************************************************************************* //tensor rtr = (C & ts ) & (C.T());
