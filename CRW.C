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

#include "CRW.H"
#include "demandDrivenData.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "KinematicParcel.H"
#include "ParticleErosion.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// this is CRW in name of CRW for making
template<class CloudType>
Foam::CRW<CloudType>::CRW
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner),
    rsPtr_(nullptr),
    ownrs_(false),
    divRsPtr_(nullptr),
    owndivRs_(false),
    gradkPtr_(nullptr),
    ownGradK_(false),
	diaPtr_(nullptr),
    owndia_(false)                           //gradkPtr_(nullptr),
                               // ownGradK_(false)
{}


template<class CloudType>
Foam::CRW<CloudType>::CRW
(
    const CRW<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm),
    rsPtr_(dm.rsPtr_),
    ownrs_(dm.ownrs_),
    divRsPtr_(dm.divRsPtr_),
    owndivRs_(dm.owndivRs_),
    gradkPtr_(dm.gradkPtr_),
    ownGradK_(dm.ownGradK_),
   // gradkPtr_(dm.gradkPtr_),
    //ownGradK_(dm.ownG#include "fvmDiv.H"radK_)
    diaPtr_(dm.diaPtr_),
    owndia_(dm.owndia_)

{
   dm.owndivRs_ = false;
    dm.ownGradK_ = false;// dm.ownGradK_ = false;
    dm.ownrs_ = false;
	dm.owndia_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRW<CloudType>::~CRW()
{
    cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CRW<CloudType>::cacheFields(const bool store)
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
	

    }

}


template<class CloudType>
Foam::vector Foam::CRW<CloudType>::update
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
    
    const scalar cps = 0.16432;

    const scalar k = this->kPtr_->primitiveField()[celli];
	const scalar nut = this->nutPtr_->primitiveField()[celli];
    const scalar epsilon =
        this->epsilonPtr_->primitiveField()[celli] + ROOTVSMALL;
	
    tensor rs = this->rsPtr_->primitiveField()[celli];
	const vector& gradk = this->gradkPtr_->primitiveField()[celli];
	const vector& divRs = this->divRsPtr_->primitiveField()[celli];
    const scalar UrelMag = mag(U - Uc - UTurb);
	const scalar dia1=this->diaPtr_->primitiveField()[celli];
    const scalar tTurbLoc =
        min(k/epsilon, cps*pow(k, 1.5)/epsilon/(UrelMag + SMALL));//UTurb = Uc/100;
   tensor ts = rs;//check if nut or mut
   ts.xx()=ts.xx()+(2.0*k/3.0);//restr calculation
	ts.yy()=ts.yy()+(2.0*k/3.0);
	ts.zz()=ts.zz()+(2.0*k/3.0);
    	vector drs = divRs;drs = drs +((2/3)*gradk) ;//divergence of re stress calculation
	//tensor rtr = (C & ts ) & (C.T()) ;
	//error passage
	vector Ur = -U + Uc + UTurb ;//formation of new coordinate system
	vector Ur1=-Ur;
	scalar ctet = Ur1.x()/sqrt((sqr(Ur1.x())+sqr(Ur1.y())));
	scalar stet= -Ur1.y()/sqrt((sqr(Ur1.x())+sqr(Ur1.y())));
	scalar cph= sqrt((sqr(Ur1.x())+sqr(Ur1.y())))/mag(Ur1);	
	scalar sph = - Ur1.z()/mag(Ur1);
	tensor C(ctet*cph,stet,ctet*sph,-stet*cph,ctet,-stet*sph,-sph,0,cph);C = C.T();
	tensor rtr = transform(C,ts) ;
	vector ri(0,0,0);
	vector tai(0,0,0);tai.x()=0.4819*1.5*rtr.xx()/epsilon;tai.y()=0.4819*1.5*rtr.yy()/epsilon;tai.z()=0.4819*1.5*rtr.zz()/epsilon;
	

	//crw start
        ri.x() = exp(-dt/(tai.x()+ROOTVSMALL) );                 //Rli calc
	ri.y() =  exp(-dt/(tai.y()+ROOTVSMALL) );
	ri.z()= exp(-dt/(tai.z()+ROOTVSMALL) );
	vector li(1,2,3);
	li.x()=(1.5*0.4819/(0.356))*((pow(rtr.xx(),1.5))/epsilon) ;
	li.y()=(1.5*0.4819/(2*0.356))*((pow(rtr.yy(),1.5))/epsilon);
	li.z()=(1.5*0.4819/(2*0.356))*((pow(rtr.zz(),1.5))/epsilon) ;  //Lei calc
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
	scalar tp =2740*(sqr(dia1))/(1.2*18*1.57e-5);scalar rerel = (dia1)*mag(Ur)/(1.57e-5); scalar coe=0;
	if (rerel<=1000)
	{coe=1+(0.15*pow(rerel,0.687));}
	else
	{coe=0.44*rerel/24 ;}
	tp =tp/coe;

     scalar sno=tp/((k/(7*epsilon))+ROOTVSMALL);//drift correction
	vector dcr = drs/(1+sno);
	vector Utr = transform(C,UTurb);       
	scalar c1 = rnd.GaussNormal<scalar>(); scalar c2 = rnd.GaussNormal<scalar>(); scalar c3 = rnd.GaussNormal<scalar>();
	vector zi(0,0,0); zi.x() = c1; zi.y() = c2; zi.z() = c3;
	UTurb = (Utr & aij ) + (zi & bij);UTurb = transform(C.T(),UTurb);
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
	aij.xx() = exp(-dt/tau.x())*exp(-dt*mag(Ur)/el.x())*(1-( (dt*mag(Ur)/(2*el.x()) ) * (1-sqr(Ur.x()/mag(Ur)))));
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




//CRW start	
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
