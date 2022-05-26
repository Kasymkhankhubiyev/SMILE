/* eint_discrt.cc
 * Internal energy generator
 * Used for generation of interal energy (or levels)
 * according to the Boltzman distribution for:
 * - freestream initiation fill
 * - freestream injection
 * - jet injection
 * - geometry reflection
 * 28.03.02         A.Kashkovsky
 * 21.06.10         A.Kashkovsky
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mytype.h"
#include "mathd.h"
#include "reader.h"
//#include "array.h"
//#include "tmatrix.h"
#include "g_gas.h"
#include "g_chem.h"
#include "col_mix.h"
#include "eint_discr.h"
#include "cells.h"
#include "particles.h"
#include "rndm.h"
#include "ele.h"

extern cgGases *Gases;        // Gases list
extern cParticles   Mol;     // Particles
extern cRndm Rn;     // random number generator
extern cCells       Cells; // Cells
extern cgMeanTv    *meanTv; // Counter mean Tv in cells

extern R8 Boltz; 

//########################################## cEIntDiscr

/*******************************************************
cEIntDiscr::cEIntDiscr
 *  constructor
 * 18.06.10       A.Kashkovsky
 *******************************************************/
cEIntDiscr::cEIntDiscr():cgEInt(),useprob(0),PR(0),PV(0)
{
  LRidx = Mol.AddMember( "LRot",'I',2*sizeof(U2), "LRot" );
  LVidx = Mol.AddMember( "LVib",'I',8*sizeof(U1), "LVib" );

}


/*******************************************************
 * cEIntDiscr::Init
 *  Initialization of discrete levels energy generation
 * 28.03.02       A.Kashkovsky
 *******************************************************/
void cEIntDiscr::Init(
 const I4 &gas,	// gas index
 R8 Trot,	// Rot. temperature
 R8 Tvib,	// Vib. Temperature
 I4 flag	// flag for use probability
)
{
 useprob = flag;
 gg=Gases->Get(gas);
 TR = Trot;
 TV = Tvib;
// calculate partition function
 gg->PartFunc( &QRot, QVib, TR, TV );
// if set useprop flag, prepere probapility
 if( useprob )
 { I4 i,l;
   // calculate dimnesion of probability for rot.
   if( gg->rotmods == 2 ) l = (gg->LevelR+1)*(gg->LevelR+2)/2;
   else                   l = gg->LevelR;
   PR = new R8[l+1]; // !! +1 - for garanty non zero only
   // calculate dimnesion of probability for vib.
   for( l=i=0;i<gg->vibmods;i++) l+=gg->LevelV[i];
   PV = new R8[l+1]; // !! +1 - for garanty non zero only
   // calculate probabilities
   gg->PartFuncProbability( PR,PV, QRot,QVib, TR,TV );
 }
}

/*******************************************************
 * cEIntDiscr::Set
 *  Generat levels of energy
 * 28.03.02       A.Kashkovsky
 *******************************************************/
void cEIntDiscr::Set(
 const I4 &pid, // particle id
 R8 Trot,       // def 0. ( for gen from flow)
 R8 Tvib
)
{
 U2 Lr[2]={0};
 U1 Lv[8]={0};
 pcgGas gas;

 if( gg ) gas = gg;
 else
 {
   I4 comp = Gases->GetPtcGasIndex(pid);
   gas = Gases->Get(comp);
 }
 if( gas->type == ATOM )
 { Lr[0]=Lr[1]=12345;
   SetE( pid, Lr, Lv );
   return;
 }

 if( useprob )
      gas->SetEnergyLevel( Lr, Lv, PR, PV);
 else
 { // calculate partition function
   if( (Trot > 1. && ABS(Trot-TR) > 1.e-4 ) ||
       (Tvib > 1. && ABS(Tvib-TV) > 1.e-4 ) )
   { TR = Trot;
     TV = Tvib;
     gas->PartFunc( &QRot, QVib, TR, TV );
   }
   gas->SetEnergyLevel( Lr, Lv, QRot, QVib, TR, TV);
 }
 SetE( pid, Lr, Lv );
}

/*******************************************************
 * cEIntDiscr::SetE
 *  Direct set internal energy levels to particle
 * 14.09.10       A.Kashkovsky
 *******************************************************/
void cEIntDiscr::SetE(
const I4 &pid,  // particle pid
U2 *Lrot,	// Rot. temperature
U1 *Lvib	// Vib. Temperature
)
{
 Mol.SetStr( pid, LRidx, (char *)Lrot );
 Mol.SetStr( pid, LVidx, (char *)Lvib );
}

/*******************************************************
 * cEIntDiscr::SetE
 *  Direct set internal energy levels to particle
 * 9.1.17  Shevyrin Al.An.
 *******************************************************/
void cEIntDiscr::SetE(
 char * p_,     // particle adress
U2 *Lrot,	// Rot. temperature
U1 *Lvib	// Vib. Temperature
)
{
  Mol.SetStr( p_,LRidx, (char *)Lrot );
  Mol.SetStr( p_, LVidx, (char *)Lvib );
}

/*******************************************************
 * cEIntDiscr::SetE
 *  Divide internal energy between modes EQUALLY
 *  (Yep, it is what we do in reactions for discret eint)
 *  like in cgGas::GetLevelsFromEnergy.
 *  Change parameters Er Ev to resulting values.
 *  25.1.17  Shevyrin Al.An.
 *******************************************************/
void cEIntDiscr::SetEInt_evenly
(const pcgGas g,// gas of ptc, just not to define
 char * p_,     // particle adress
 R8 &Er,	// total Erot 
 R8 &Ev	        // total Evib
 )
{
  U2 Lr[2]={0};	// Rotational energy levels
  U1 Lv[8]={0};	// Vibrational energy levels
  g->GetLevelsFromEnergy(Er, Ev, Lr, Lv, &Er, &Ev); // problem?
  Mol.SetStr( p_, LRidx, (char *)Lr );
  Mol.SetStr( p_, LVidx, (char *)Lv );
} // end cEIntDiscr::SetEInt_evenly



/*******************************************************
 * cEIntDiscr::GetE
 *  Direct extract internal energy levels
 *  from particle
 * 14.09.10       A.Kashkovsky
 *******************************************************/
void cEIntDiscr::GetE(
const I4 &pid,  // particle pid
U2 *Lrot,	// Rot. temperature
U1 *Lvib	// Vib. Temperature
)
{

 Mol.GetStr( pid, LRidx, (char *)Lrot );
 Mol.GetStr( pid, LVidx, (char *)Lvib );

}

/*******************************************************
 * cEIntDiscr::GetE
 *  Extract and calculate internal energy 
 *  from particle
 * 14.09.10       A.Kashkovsky
 * replaced by direct get
 * 29.12.16  Shevyrin Al.An.
 *******************************************************/
void cEIntDiscr::GetE(
const I4 &pid, // particle pid
 R8 *ER,       // Rot. temperature
 R8 *EV       // Vib. Temperature
)
{
  char* p = Mol.pMolGet(pid);
  GetE_ptr(p,ER,EV);
}

/*******************************************************
 * cEIntCont::GetE
 *  Direct extract internal energy from particle
 * 29.12.16  Shevyrin Al.An.
 *******************************************************/
void cEIntDiscr::GetE_ptr(
 char * p_,     // particle adress
 R8 *Erot,	// Rot. temperature
 R8 *Evib	// Vib. Temperature
)
{
 U2 Lrot[2]; // Rot. energy levels
 U1 Lvib[8]; // Vib. energy levels

 Mol.GetStr( p_, LRidx, (char *)Lrot );
 Mol.GetStr( p_, LVidx, (char *)Lvib );

 I4 comp; 
 if (Gases->Gidx <0 ) comp=0; 
 else Mol.Get( p_, Gases->Gidx, &comp ); // Set gas
 pcgGas g = Gases->Get( comp );     // get gas ptr

 (*Erot)=(*Evib)=0.;
 if( g->type != ATOM ){
   I4 j;
   // rot. energy
   for(j=0;j<g->rotmods;++j) (*Erot)+=g->ER[j] [Lrot[j] ];
   // vib. energy
   for(j=0;j<g->vibmods;++j) (*Evib)+=g->EV[j] [Lvib[j] ];
 }

}

/*******************************************************
 * cEIntDiscr::GetEMod
 *  Extract and calculate internal energy 
 *  from particle with output by vibration modes
 * 14.09.10       A.Kashkovsky
 *******************************************************/
void cEIntDiscr::GetEMod(
const I4 &pid, // particle pid
 R8 *ER,       // Rot. temperature
 R8 *EV       // Vib. Temperature [8]!!!!!!
)
{
 U2 Lrot[2]; // Rot. energy levels
 U1 Lvib[8]; // Vib. energy levels

 Mol.GetStr( pid, LRidx, (char *)Lrot );
 Mol.GetStr( pid, LVidx, (char *)Lvib );
 I4 comp = Gases->GetPtcGasIndex(pid);

 pcgGas g = Gases->Get(comp);
 (*ER)=EV[0]=EV[1]=EV[2]=EV[3]=EV[4]=EV[5]=EV[6]=EV[7]=0.;
 if( g->type != ATOM )
 { I4 j;
  // rot. energy
   for( j=0;j<g->rotmods;j++) (*ER)+=g->ER[j] [Lrot[j] ];
  // vib. energy
   for( j=0;j<g->vibmods;j++) EV[j] =g->EV[j] [Lvib[j] ];
 }

}

/*******************************************************
 * cEIntCont::GetEMod
 *  multiply eint by factor.
 *  use modes in molecules as base values
 * 24.1.17. Shevyrin Al.An.
 * do for ptc[r/w]
 * 12.2.19 Shevyrin Al.An.
 *******************************************************/
void cEIntDiscr::ChangeEintBy_m
(char* ptcr,        // particle read
 char* ptcw,        // particle write
 const pcgGas g,
 const R8 &factor, // factor of change
 R8 &resER,        // returing result of operations
 R8 &resEV         // returing result of operations
 )
{
 U2 Lrot[2]; // Rot. energy levels
 U1 Lvib[8]; // Vib. energy levels
 Mol.GetStr( ptcr, LRidx, (char *)Lrot );
 Mol.GetStr( ptcr, LVidx, (char *)Lvib );

// E rot
  resER=0.;
  for(I4 i=0;i<g->rotmods;i++)
  { 
    R8 E = factor * g->ER[i][ Lrot[i] ];
    R8 LmaxR = 0.5*(-1.+sqrt( 1.+4.*E/(Boltz*g->TetaR[i]) ));
    I4 Lmax=(I4)(LmaxR+0.001);
    if( g->ER[i][Lmax] > E ) --Lmax;
    if( Lmax < 0 ) Lmax = 0;
    resER+=g->ER[i][Lmax];
    Lrot[i] = Lmax;
  }

// E vib
  resEV=0.;

  for(I4 i=0;i<g->vibmods;i++)
  { R8 E = factor * g->EV[i][ Lvib[i] ] + g->EV0[i];
    R8 LmaxR = (g->EtaV[i] > 0.0001 ?
    /*	
     Ev0 =  0.5*Boltz*mix->gL->TetaV[i]*
       ( 1.-mix->gL->EtaV[i]*Boltz*mix->gL->TetaV[i]/
               (8.*mix->gL->EDiss) );
*/
		 2.*g->EDiss/(Boltz*g->TetaV[i])*
		( 1.-sqrt(1.-g->EtaV[i]*E/g->EDiss) ) /
		g->EtaV[i]
		-0.49999 // it isn't 0.5, because 13.99999999 == 14
		: E/(Boltz*g->TetaV[i]) - 0.49999);
    I4 Lmax=(I4)LmaxR;
    if( g->EV[i][Lmax] > E ) Lmax--;
    if( Lmax < 0 ) Lmax = 0;
    resEV+=g->EV[i][Lmax];
    Lvib[i] = Lmax;
  }

  Mol.SetStr( ptcw, LRidx,(char *)Lrot );
  Mol.SetStr( ptcw, LVidx,(char *)Lvib );
}


/********************************************************
 * cEIntDiscr::GetVibDOFfromEnergy
 *  Return Vib. degree of freedom, calculated from energy
 *  For discrete model
 * 09.04.02			Kashkovsky A.V.
 * 06.04.06			Kashkovsky A.V.
 * 14.09.10			Kashkovsky A.V.
 ********************************************************/
R8 cEIntDiscr::GetVibDOFfromEnergy(
 const I4 &pid     // particle index
)
{
 I4 comp = Gases->GetPtcGasIndex(pid);
 pcgGas g = Gases->Get(comp);
 if( g->type == ATOM ) return(0.);
 R8 e;
 R8 dof = 0.;
/***************
The follow formules (for Simple Harmonic Oscilator) is used:
 $$
\zeta_v = \sum_m^{vibmodes} \frac{2\Theta_m/T}{exp(\Theta_m/T)-1}
   \cdot dgen_m
$$
$$ E = \frac{\zeta}{2}kT $$
and easy follow:
$$
 \zeta_v = \sum_m^{vibmodes} \frac{2E}{k\Theta_m}\cdot
 log(\frac{k\Theta_m}{E}+1)  \cdot dgen_m
$$

--+ THIS FORMYLA HAS TO BE REMOVE!! becouse we don't use SHO
--+ It has to be replase by calculation dof from Energy level
--+ by interpolation
--+ 04.02.03

****************/
 int mv(-1);

 U1 Lvib[8]; // Vib. energy levels

 Mol.GetStr( pid, LVidx, (char *)Lvib );

 if( mv < 0 ) // all modes
 { I4 i;
   for( i=0;i<g->vibmods;i++ )
   {
     e = g->EV[i][ Lvib[i] ]/(Boltz*g->TetaV[i]);
     if( e < 1.e-10 ) continue;
     dof+= 2.*e*log( 1.+1./e )*g->degevib[i];
   }
 }
 else
 {
   e = g->EV[mv][ Lvib[mv] ]/(Boltz*g->TetaV[mv]);
   if( e < 1.e-10 ) return(0.);
   dof+= 2.*e*log( 1.+1./e )*g->degevib[mv];
 }
 return(dof);
}

/********************************************************
 * cEIntDiscr::GetVibDOFfromEnergy
 *  Return Vib. degree of freedom, calculated from energy
 *  For discret model
 * 02.04.14       Shevyrin Al.An.
 ********************************************************/
R8 cEIntDiscr::GetVibDOFfromEnergy(
 const pcgGas g, // component of gas
 const R8 &Evib  // vibration energy
)
{
 if( g->type == ATOM ) return(0.);
 R8 dof = 0.;
 for( I4 i=0;i<g->vibmods;++i ) {
   const R8 e = Evib/(Boltz*g->TetaV[i]);
   if( e < 1.e-10 ) continue;
   dof+= 2.*e*log( 1.+1./e )*g->degevib[i];
 }
 return(dof);
}



/********************************************************
 * cEIntDiscr::GetVibDOF_TCE
 *  Return Vib. degree of freedom, calculated from temperature
 * 22.5.12          Shevyrin Al.An.
 * use exp() expression instead of log()
 * 3.5.14           Shevyrin Al.An.
 ********************************************************/
R8 cEIntDiscr::GetVibDOF_TCE(const pcgGas g,const R8 &T){
  if( g->type == ATOM ) return(0.);
  R8 dof = 0.;
  for(I4 i=0;i<g->vibmods;++i){
    R8 e = g->TetaV[i]/T;
    if( e < 40. )
      dof+= 2.*e/(exp(e)-1.)*g->degevib[i];
  }
  return dof;
};

/*******************************************************
 * cEIntDiscr::SetEnergy
 *  set int.e. to molecule,  does not define gas
 *  3.4.14. Shevyrin Al.An.
 *******************************************************/
void cEIntDiscr::SetEnergy(const pcgGas g,const I4 &pid, R8 *Erot, R8 *Evib){
  U2 Lr[2]={0};	// Rotational energy levels
  U1 Lv[8]={0};	// Vibrational energy levels
  R8 ErL(*Erot), EvL(*Evib);// for input
  g->GetLevelsFromEnergy(ErL, EvL, Lr, Lv, Erot, Evib );
  SetE(pid, Lr, Lv); // only 25.1.17  Shevyrin Al.An. added
}; // end cEIntDiscr::SetEnergy(const pcgGas,const I4&,R8*,R8*)


/*******************************************************
 * cEIntDiscr::EintExch
 *  Exchange of internal energy
 * 14.09.10       A.Kashkovsky
 *******************************************************/
void cEIntDiscr::EintExch(
 cgCol *col // collision class
)
{
/*================ from cgCol class (g_col.h):

  I4 pL;	// Index of ptc. in particle container
  I4 pM;
  pcgGas gL;	// Gas description
  pcgGas gM;
  V8 vL;       // velocity of particle
  V8 vM;
  V8 Vrel;	// Relative velocity
  R8 Vmod;	// Module  of Relative velocity
=============*/
 R8 Erel;	// Relative energy
 R8 Rnd,Dust,F,W;
 U2 Lr[2];	// Rotational energy levels
 U1 Lv[8];	// Vibrational energy levels
 I4 m;

 if( col->gL->type == ATOM && col->gM->type == ATOM ) return;

// check if there is internal energy transfer according to LB scheme
 if( Rn.dm() > ProbSum.Data(col->igL, col->igM) ) return;

 Erel = col->gL->c_massr[col->igM]*col->Vmod2*0.5;

// check with first molecule 
 Rnd = Rn.dm();
 Dust = col->gL->rotmods*ProbRT.Data( col->igL, col->igM)/
                    ProbSum.Data(col->igL, col->igM);
 if( Rnd < Dust )
 { // RT transfer for L molecule
   Mol.GetStr( col->ptc1L, LRidx, (char *)Lr );
   //   GetE_ptr(col->ptc1L, Lr, Lv );
   ColRot(col->gL,col->gM, Lr, Erel );
   Mol.SetStr( col->ptc2L, LRidx, (char *)Lr );
   //SetE_ptr(col->ptc2L, Lr, Lv );
   goto END;
 }
 Rnd-=Dust;
 for( m=0;m<col->gL->vibmods;m++)
 { Dust = ProbVT[m].Data(col->gL->id, col->gM->id)/
          ProbSum.Data(col->gL->id, col->gM->id);
   if( Rnd < Dust )
   { // VT transfer for L molecule
     Mol.GetStr( col->ptc1L, LVidx, (char *)Lv );
     //      GetE_ptr(col->ptc1L, Lr, Lv );
     ColVib(col->gL,col->gM, m,Lv, Erel );
     Mol.SetStr( col->ptc2L, LVidx, (char *)Lv );
     //     SetE_ptr(col->ptc2L, Lr, Lv );
     goto END;
   }
   Rnd-=Dust;
 } // end for( m=0;m<col->gL->vibmods;m++)

// check with second molecule 
 Dust = col->gM->rotmods*ProbRT.Data(col->igM, col->igL)/
                    ProbSum.Data(col->igM, col->igL);
 if( Rnd < Dust )
 { // RT transfer for M molecule
   Mol.GetStr( col->ptc1M, LRidx, (char *)Lr );
   //   GetE_ptr(col->ptc1M, Lr, Lv );
   ColRot(col->gM,col->gL, Lr, Erel );
   Mol.SetStr( col->ptc2M, LRidx, (char *)Lr );
   //   SetE_ptr(col->ptc2M, Lr, Lv );
   goto END;
 }
 Rnd-=Dust;
 for( m=0;m<col->gM->vibmods;m++)
 { Dust = ProbVT[m].Data(col->igM, col->igL)/
          ProbSum.Data(col->igM, col->igL);
   if( Rnd < Dust )
   { // VT transfer for M molecule
     Mol.GetStr( col->ptc1M, LVidx, (char *)Lv );
     //     GetE_ptr(col->ptc1M, Lr, Lv );
     ColVib(col->gM,col->gL, m,Lv, Erel );
     Mol.SetStr( col->ptc2M, LVidx, (char *)Lv );
     //     SetE_ptr(col->ptc2M, Lr, Lv );
     goto END;
   }
   Rnd-=Dust;
 } // end for( m=0;m<col->gM->vibmods;m++)

END:
// change relative velocity 
 W = sqrt(2.* Erel/ col->gL->c_massr[col->igM]);
 F = W/col->Vmod;
 col->Vrel*=F;
 col->Vmod =W;


 
}

/*******************************************************
 * cEIntDiscr::ColRot
 *  Rotational energy exchange for Discrete model.
 * 4.4.2			Kashkovsky A.V.
 * 14.09.10       A.Kashkovsky
 *******************************************************/
void cEIntDiscr::ColRot(
 pcgGas gA,	// gas property for A particle
 pcgGas gB,	// gas property for B particle
 U2 *Lr,	// Levels of rotational energy
 R8 &Erel	// Relative energy
)
{
 R8 Ecol;     // energy of collision
 R8 Jmax;     // max. Rotational level
 I4 II;       // level to be choised
 I4 JJ;       // level where max. of distribution functin
 R8 Fmax,F,D; // Fmax=f(JJ)
 R8 Err;      // energy to pass to rot energy. 
 R8 ff;       // temporary variable

 switch( gA->type )
 {
 case ATOM: break;
 case DIATOMIC:
 case LINEAR:
      Ecol = Erel + gA->ER[ 0 ][ Lr[0] ];
      Err = Ecol < gA->EDiss ?  Ecol : gA->EDiss;
      /* 16.11.10 sasa
        Find  max. possible level for energy Err
         
         For rigid rotator, which is used here
         $$
         J_{max}=\left[0.5\cdot\left(-1+\sqrt{1+{4E'\over
         k\theta_r}}\right)\right]
         $$
      */
      Jmax = ( 0.5*(-1.+sqrt( 1.+4.*Err/(Boltz*gA->TetaR[0])
                                   ) ) );
      // level has to be choice uniformly
      // round Jmax to hight integer, for make 
      // probability of choice last level the same,
      // as for all others level
      Jmax = (R8)((int)Jmax+1) - 1.e-9;
      if( Jmax < 1.0 )
      { Lr[0] = 0;
        Erel = Ecol;
        return;
      }
      /* 19.02.2016 sasa
        ER[Jmax] was a little bit higher Ecol => crash
        Check this.
      */ 
      if( gA->ER[0][(int)Jmax] > Err ) Jmax -= 1.0;
      /* 16.11.10 sasa
         find a level, for which the distribution function
         is a max.
         $$
         J^\star=0.5\cdot
            \left(-1+\sqrt{{4E'\over k\theta_r}+1\over
                    3-2\alpha}\right)
         $$
      */
      JJ = (int)( 0.5*(-1.+sqrt( (1.+4.*Err/(Boltz*gA->TetaR[0]))/
                                  (3.-2.*gA->c_alpha[gB->id])
                               ) ) );
      if( JJ > gA->LevelR-1 ) JJ = gA->LevelR-1;
      //  distrb. function for level JJ
      Fmax= (2.* JJ   +1.)*
             pow( Ecol-gA->ER[0][JJ],   1.-gA->c_alpha[gB->id] );

      // the same, for JJ+1 (one of them is maximum
      if( JJ < gA->LevelR-1 && Ecol > gA->ER[0][JJ+1])
      { F = (2.*(JJ+1)+1.)*
             pow( Ecol-gA->ER[0][JJ+1], 1.-gA->c_alpha[gB->id] );
        if( F > Fmax ) Fmax = F;
      }
      // choice a level
      while(1)
      { II = (int)( Jmax*Rn.dm() );
        F = (2.* II   +1.)*
             pow(Ecol-gA->ER[0][II], 1.-gA->c_alpha[gB->id])/
            Fmax;
        if( Rn.dm() < F ) break; 
      }
      Lr[0] = II;
      Erel = Ecol - gA->ER[ 0 ][ II ];
      break;
 case SPHERICAL:
      Ecol = Erel + gA->ER[ 0 ][ Lr[0] ];
      Err = Ecol < gA->EDiss ?  Ecol : gA->EDiss;
      /* 16.11.10 sasa
        Find  max. possible level for energy Err
        the same as for LINEAR 
      */
      Jmax = ( 0.5*(-1.+sqrt( 1.+4.*Err/(Boltz*gA->TetaR[0])
                                   ) ) );
      // level has to be choice uniformly
      // round Jmax to hight integer, for make 
      // probability of choice last level the same,
      // as for all others level
      Jmax = (R8)((int)Jmax+1) - 1.e-9;
      if( Jmax < 1.0 )
      { Lr[0] = 0;
        Erel = Ecol;
        return;
      }
      JJ = (int)( 0.5*(-1.+sqrt( (1.+4.*Err/(Boltz*gA->TetaR[0]))/
                                  (2.-gA->c_alpha[gB->id])
                               ) ) );
      if( JJ > gA->LevelR-1 ) JJ = gA->LevelR-1;
      //  distrb. function for level JJ
      ff = (2.* JJ   +1.);
      Fmax= ff*ff*
             pow( Ecol-gA->ER[0][JJ],   1.-gA->c_alpha[gB->id] );

      // the same, for JJ+1 (one of them is maximum
      if( JJ < gA->LevelR-1 && Ecol > gA->ER[0][JJ+1])
      { ff = (2.*(JJ+1)+1.);
        F = ff*ff*
             pow( Ecol-gA->ER[0][JJ+1], 1.-gA->c_alpha[gB->id] );
        if( F > Fmax ) Fmax = F;
      }
      // choice a level
      while(2)
      { II = (int)( Jmax*Rn.dm() );
        ff = (2.*II+1.);
        F =  ff*ff*
               pow(Ecol-gA->ER[0][II],    1.-gA->c_alpha[gB->id])/
                                                             Fmax;
        if( Rn.dm() < F ) break; 
      }
      Lr[0] = II;
      Erel = Ecol - gA->ER[ 0 ][ II ];
      break;
 case GENERAL:
      if( Rn.dm() < 0.5 )
      { // change  first mode
        Ecol = Erel + gA->ER[ 0 ][ Lr[0] ];
        Err = Ecol < gA->EDiss ?  Ecol : gA->EDiss;
        /* 16.11.10 sasa
          Find  max. possible level for energy Err
          the same as for LINEAR 
        */
        Jmax = ( 0.5*(-1.+sqrt( 1.+4.*Err/(Boltz*gA->TetaR[0])
                                    ) ) );
        Jmax = (R8)((int)Jmax+1) - 1.e-9;
        if( Jmax < 1.0 )
        { Lr[0] = 0;
          Erel = Ecol;
          return;
        }
        
        JJ=(int)( 0.5*(-1.+sqrt( (1.+4.*Err/(Boltz*gA->TetaR[0]))/
                                  (3.-2.*gA->c_alpha[gB->id])
                               ) ) );
//        if( JJ > gA->LevelR-1 ) JJ = gA->LevelR-1;
        if( JJ > Jmax ) JJ = (I4) Jmax;
        //  distrb. function for level JJ
        Fmax= (2.* JJ   +1.)*
               pow( Ecol-gA->ER[0][JJ],   1.-gA->c_alpha[gB->id] );

        // the same, for JJ+1 (one of them is maximum
        if( JJ < gA->LevelR-1 )
        { F = (2.*(JJ+1)+1.)*
               pow( Ecol-gA->ER[0][JJ], 1.-gA->c_alpha[gB->id] );
          if( F > Fmax ) Fmax = F;
        }
        // choice a level

        while(3)
        { II = (int)( ((R8)(Jmax+1-Lr[1]))*Rn.dm() )+Lr[1];
          if( II > Jmax ) II = (I4) Jmax;
          F=   (2.* II   +1.)*
             pow(Ecol-gA->ER[0][II],    1.-gA->c_alpha[gB->id])/Fmax;
          if( Rn.dm() < F ) break; 
        }
        Lr[0] = II;
        Erel = Ecol - gA->ER[ 0 ][ II ];
      }
      else // if( Rndm() < 0.5 )
      { // change of second mode
        Ecol = Erel + gA->ER[ 1 ][ Lr[1] ];
        Fmax = pow( Ecol, 1.-gA->c_alpha[gB->id]);
        F = 2.*pow( Ecol-gA->ER[1][0], 1.-gA->c_alpha[gB->id]);
        if( F > Fmax ) Fmax = F;
        while(4)
        { II = (int)( ((R8)Lr[0]+1.)*Rn.dm() );
          if( II ) D = 2.;
          else     D = 1.;
          F  = pow( D*(Ecol-gA->ER[1][II]), 1.-gA->c_alpha[gB->id])/
                                                            Fmax;
          if( Rn.dm() < F ) break; 
        }
        Lr[1] = II;
        Erel = Ecol - gA->ER[ 1 ][ II ];
      }
      break;
 }
}

/********************************************************
 * cEIntDiscr::ColVib
 *  Vibrational energy exchange for Discrete model.
 * 5.4.2			Kashkovsky A.V.
 * 14.09.10 A.Kashkovsky 
 * 17.11.10 A.Kashkovsky Formulas according theory.pdf
 * 19.05.20 A.Latyshev for polyatomic molecules total energy is less than Ed
 ********************************************************/
void cEIntDiscr::ColVib(
 pcgGas gA,	// gas property for A particle
 pcgGas gB,	// gas property for B particle
 I4  mode,	// mode number
 U1 *Lv,	// Levels of vibrational energy
 R8 &Erel	// Relative energy
)
{
 if( gA->type == ATOM )
 { Lv[mode] = 0;
   return;
 }

 R8 Ecol; 	// energy of collision
 R8 Vmax;	// max. Vib. level
 I4 II;         // level index

/* 16.11.10 sasa
 Choice of Vmax wasre-made according "theory.pdf":

For anharmonic oscillator 
$$
V_{max}=\left[{2E_d\over k\theta_v}\cdot {1-\sqrt{1-\eta E'/E_d} \over
\eta}\right], \;\; E'=min(E_c, E_d) 
$$

For harmonic oscillator ($\eta=0$),
$$
V_{max} = \left[{E_c \over k\theta_v}\right],
$$

$$
F(V)/F_{max} = (E_c-E_v(V))^{1-\alpha}/E_c^{1-\alpha} =
 (1-E_v(V)/E_c)^{1-\alpha}
$$

$$
E_v(V)=k\theta_v(V+{1\over 2})
$$


*/

// total energy of vibration
R8 Eallv = 0;
for (int m=0; m<gA->vibmods; m++)
{
  Eallv += gA->EV[m][ Lv[m] ];
}

// total energy for split over Trans. and Vib
 Ecol = Erel + gA->EV[mode][ Lv[mode] ];
// energy, which has to be fit in vib. mode
 R8 Evv = (Eallv+Erel) < gA->EDiss ? Ecol : (gA->EDiss-(Eallv-gA->EV[mode][ Lv[mode] ]));

//We ban anharmonicity for polyatomic molecules
 if(( gA->EtaV[mode] > 0.0001 ) & ( gA->type == DIATOMIC ))
 { // anharmonic oscillator
   Vmax = 2.*gA->EDiss/(Boltz*gA->TetaV[mode])*
           ( 1.-sqrt(1.-gA->EtaV[mode]*Evv/gA->EDiss)
           ) /  gA->EtaV[mode];
 }
 else // harmonic oscillator
 { Vmax = Evv/(Boltz*gA->TetaV[mode]);
 }
// level has to be choice uniformly
// round Vmax to hight integer, for make 
// probability of choice last level the same,
// as for all others level
 Vmax = (R8)((int)Vmax+1) - 1.e-9;
 if( Vmax > gA->LevelV[mode] ) Vmax = gA->LevelV[mode];

// collision energy less then energy of first level
// no enegry in vibration
 if( Vmax < 1.0 )
 { Lv[mode] = 0;
   Erel = Ecol;
   return;
 }

// choice a level 
 R8 F; // distribution function result for level II

// R8 epow,dm,dmm;
 while(1)
 { II = (int)( Vmax*Rn.dm() );
/*
   epow = gA->EV[mode][II]/Ecol;
   dm  = gA->DegenMult(mode,II);
   dmm = gA->DegenMult(mode,(I4)(Vmax+1.));
   F  = pow( 1.- epow,1.-gA->c_alpha[gB->id] )*dm/dmm;
*/
   F  = pow( 1.- gA->EV[mode][II]/Ecol,
             1.- gA->c_alpha[gB->id] 
           )*gA->DegenMult(mode,II)/
             gA->DegenMult(mode,(I4)(Vmax+1.));

   if( Rn.dm() < F ) break;    
 }

 Lv[mode] = II;
 Erel = Ecol - gA->EV[mode][ II ]; 

}
// #################### 


/********************************************************
 * cEIntDiscr::ZrZvUpdate
 *  Internal Energy exchange, 
 * Recalculate  ZrZv Variable dependent from T
 * 18.02.11			Kashkovsky A.V.
 ********************************************************/
void cEIntDiscr::ZrZvUpdate(
 R8 T		// temperature
)
{
  
/*
  $Z_r$ and $Z_v$ in DSMC is differ from those used in continum  
$Z_r^C$ and $Z_v^C$:

$$
Z_r = \frac{\zeta_t}{\zeta_{tr}}Z_r^C
$$
$$
Z_v = \frac{\zeta_t}{\zeta_{trv}}Z_v^C
$$

$$
\zeta_t = 4-2\alpha ; \;\;\;\;\;
\zeta_{tr} = \zeta_t + \zeta_r ; \;\;\;\;\;
\zeta_{trv} = \zeta_t + \zeta_r + \zeta_v
$$

For constant Zr Zv it's enought
*/
 int i,j,m;
 pcgGas g;
 R8 zeta_t, zeta_tr, zeta_trv, dofv,ff;
   
 if( !Gases->VarZrZv )  // T-undependent ZrZv
 { for( i=0;i<Gases->Num();++i )
   { g = Gases->Get(i);
     if( g->type == ATOM ) continue;
    
     for( j=0;j<Gases->Num();++j )
     { 
       // d.o.f
       zeta_t  = 4.-2.*g->c_alpha[j];
       zeta_tr = zeta_t + g->dofr;
       Zr.Data(i,j) = g->c_Zr[j]*zeta_t/zeta_tr;

       for(m=0;m<Gases->vibmods;++m)
       { ff = g->TetaV[m]/T;
         if( ff < 25. ) dofv = 2.*ff/(exp(ff)-1.);
         else           dofv = 0.;
         zeta_trv = zeta_tr + dofv;
         Zv[m].Data(i,j) = g->c_Zv[j]*zeta_t/zeta_trv;
       }
     }
   }
 } // end if( !Gases->VarZrZv )
 else // T-dependent 
 {
/**
Temperature dependent:

$$
Z_r^C(T) = \frac{Z_{r,\infty}}
{1 + \frac{\pi^{3/2}}{2} 
\left(\frac{T^\star}{T}\right)^{1/2} + 
\left(\frac{\pi^2}{4} + \pi\right)\frac{T^\star}{T}}
$$
$T^\star$ and $Z_{r,\infty}$ - gas dependent (input)

$$
Z_v(T) = \frac{\zeta_t}{\zeta_t+A} \frac{\tau_v}{\tau_c}
$$

$$
A = \frac{exp(\Theta_v/T)}{2}\left(\frac{2 \Theta_v/T}{exp(\Theta_v/T) -1 }\right)^2
$$
$\Theta_v$ - for component, for which $Z_v$ is calculated

$$
\frac{\tau_v}{\tau_c} = C_{VT1}T^{-0.5-\alpha}(\exp(C_{VTA} T^{-1/3})\cdot C_{VTB} + C_{VT2} \sqrt{T}) 
$$
$$
C_{VTA} = A^\star \;\;\;\;
C_{VTB} = \frac{\exp(B^\star) 101325}{  k}
$$
$$
C_{VT1} = 2d_{ref}^2 
\sqrt{\frac{\pi 2k}{m_r}} T_{ref}^\alpha
\;\;\;\;
C_{VT2} = \sqrt{\frac{\pi m_r}{4 k}}\cdot 10^{20}
$$


*/
   R8 Zrc,TvTc,AA,ThetaT;
   for( i=0;i<Gases->Num();++i )
   { g = Gases->Get(i);
     if( g->type == ATOM ) continue;

     Zrc = g->Zroo/(1.+g->ConRT1/sqrt(T)+g->ConRT2/T);

     for( j=0;j<Gases->Num();++j )
     { 
       // d.o.f
       zeta_t  = 4.-2.*g->c_alpha[j];
       zeta_tr = zeta_t + g->dofr;
 
       Zr.Data(i,j) = Zrc*zeta_t/zeta_tr;

       for(m=0;m<Gases->vibmods;++m)
       {
         if( m >= g->vibmods  ) // sasa 12.04.2012, shevr 5.9.2014
         { Zv[m].Data(i,j) = 0.;
           continue;
         }
         TvTc = Gases->ConVT1[m].Data(i,j)* 
                  pow( T, -0.5-g->c_alpha[j] ) *
             (
              exp(Gases->ConVTA[m].Data(i,j)*
                      pow( T, -1./3.) 
             ) *Gases->ConVTB[m].Data(i,j) +
              Gases->ConVT2[m].Data(i,j) *sqrt(T)
             );
         ThetaT = g->TetaV[m]/T;
         if( ThetaT < 25. ) {ff=exp(ThetaT);dofv = 2.*ThetaT/(ff-1.);}
         else               {ff=dofv = 0.;}
         AA = ff/2.*dofv*dofv;
         Zv[m].Data(i,j) = zeta_t/(zeta_t+AA)*TvTc;
       } // end for(m=0;m<Gases->vibmods;++m)
     }// end for( j=0;j<Gases->Num();++j )
   } // end for( i=0;i<Gases->Num();++i )
 } // end T-dependent 

}


// probability of Int.enegry exchange

/********************************************************
 * cEIntDiscr::ProbIntEExchage
 *  Calculation Cross Collision Probability for
 *  Internal Energy exchange, 
 * 22.02.11			Kashkovsky A.V.
 ********************************************************/
void cEIntDiscr::ProbIntEExchage()
{
 I4 i,j,m;
 pcgGas g,p;
 
// **** calcultion Rot. And Vib. probability

// set 0 values
 ProbRT.Set(0.);
 for(m=0;m<Gases->vibmods;++m) ProbVT[m].Set(0.);
 ProbSum.Set(0.);

 for( i=0;i<Gases->Num();++i )
 { g = Gases->Get(i);
   if( g->type == ATOM ) continue;

   // Rotation
   for( j=0;j<Gases->Num();++j )
   { if( Zr.Data(i,j) < 1.e-6 ) continue;
     ProbRT.Data(i,j) = 1./Zr.Data(i,j);
   }
     
   // Vibrational
   for(m=0;m<g->vibmods;++m)
   { for( j=0;j<Gases->Num();j++ )
     { if( Zv[m].Data(i,j) < 1.e-6 ) continue;
       ProbVT[m].Data(i,j) = 1./Zv[m].Data(i,j);
     }
   } // end for(m=0;m<g->vibmods;++m)
 } // end for( i=0;i<Gases->Num();++i )


// **** calculation Sum probability

 for( i=0;i<Gases->Num();i++)
 {p =  Gases->Get(i);
   for( j=0;j<Gases->Num();j++)
   { g =  Gases->Get(j);
   
     ProbSum.Data(i,j)+=p->rotmods*ProbRT.Data(i,j);
     for(m=0;m<p->vibmods;m++)
        ProbSum.Data(i,j)+= ProbVT[m].Data(i,j);
   
     ProbSum.Data(i,j)+=g->rotmods*ProbRT.Data(j,i);
     for(m=0;m<g->vibmods;m++)
        ProbSum.Data(i,j)+= ProbVT[m].Data(j,i);
   
   } // end for( j=0;j<=i;j++)
 } // end for( i=0;i<Num();i++)

}



// #################### Chemical reactions

/*******************************************************
 * cEIntDiscr::RealizeRecombination
 *  Realisation of recombination reaction
 * reaction of type A + B + M -> AB + M
 * A,B -- unstructured.
 * M -- is orbitrary particle, stored in (*InnerProb)
 * Internal energy of M does not change
 * 26.12.11 Shevyrin Al.An.
 *******************************************************/
void cEIntDiscr::RealizeRecombination(
 cMix *mix,	// mixture collision 
 cgReacRecombination* Reac // reaction 
)
{
/***********
 general idea
  molecules A & B stored in mix. 
  Third body M stored in Reac->InnerProb.
 1) calculate total energy
 2) find d.o.f., ErL, EvL .
 3) change level of int. energy of L by factor for each mode
 4) change translation energy correspondent a change of int. eng
 5) make collision
 6) delete source particles
 7) Set parameters of recombination product AB (stored in mix L)
 8) Remove B molecule
 9) Set velocity of third body (Molecule M)
************/



  /*  
      some helpfull maths:
      $$m_1 c_1^2+m_2 c_2^2 + m_3 c_3^2= (m_1+m_2)c_m^2 + m_r c_r^2 + m_3 c_3^2$$
      now quasi-particle has mass $(m_1 + m_2)$ and velocity $c_m$
      $$c_{r3}=c_m-c_3$$
      $$c_{m3}={ (m_1+m_2) c_m + m_3 c_3 \over m_1+m_2+m_3 } 
      ={  m_1 c_1 + m_2 c_2 + m_3 c_3 \over m_1+m_2+m_3} $$
      $$m_{r3} = { m_r m_3 \over m_r+m_3 }$$
  */
  const R8 Etotal = Reac->InnerProb->Ec3+Reac->Qr;

  const pcgGas P1 = Reac->P[0]; // gas of Product 1 (AB)
  const pcgGas P2 = Reac->InnerProb->gX; // gas of Product 2 (M)

// find d.o.f.
  const R8 doft = 4.-2.*P1->c_alpha[P2->id]; // trans. d.o.f
  const R8 Ttr = Etotal/(doft*Boltz);	// temp. from trans. energy
  R8 dofv1=0.;
  for(I4 i=0;i<P1->vibmods;i++)
  { const R8 TetaTr = P1->TetaV[i]/Ttr;
    if( TetaTr < 40. )
    { dofv1+=2.*TetaTr/
           (exp(TetaTr)-1.);
    }
  }
  const R8 dof = doft+P1->dofr+dofv1;
  R8 ErL = P1->dofr/dof * Etotal;
  R8 EvL = dofv1   /dof * Etotal;
  R8 Er1, Ev1;
  U2 Lr1[2]={0};	// Rotational energy levels
  U1 Lv1[8]={0};	// Vibrational energy levels
  P1->GetLevelsFromEnergy(ErL, EvL, Lr1, Lv1, &Er1, &Ev1 );
  const R8 Etr = Etotal - Er1 - Ev1;
  mix->Vmod2 = 2.*Etr/P1->c_massr[P2->id];
  mix->Vmod = sqrt(mix->Vmod2);
// simulate relative velocity for new particles
  const R8 F = mix->Vmod ,
    C = 1.-2.*Rn.dm(),
    S = F*sqrt(1.-C*C),
    Theta = 2.*M_PI*Rn.dm();
  const V8 Vr(F*C,S*cos(Theta),S*sin(Theta));

  if( mix->oL > mix->oM )
  { mix->DelPtcInCurCell( mix->oL, mix->gL->id);
    mix->DelPtcInCurCell( mix->oM, mix->gM->id);
  }
  else
  { mix->DelPtcInCurCell( mix->oM, mix->gM->id);
    mix->DelPtcInCurCell( mix->oL, mix->gL->id);
  }

  // Set parameters of recombination product (Molecule AB)
  mix->gL = P1;  //  convert A to AB molecule 
  Gases->SetPtcGasIndex( mix->pL,  &mix->gL->id);
  mix->vL = Reac->InnerProb->Vm3 + Vr*P2->c_mass[P1->id];
  Mol.Set( mix->pL, mix->Vidx, mix->vL );
  // set internal energy
  SetE( mix->pL, Lr1, Lv1 );
  // re-index particle in the cell
  mix->AddPtcInCurCell( mix->pL,  mix->gL->id);

  // Remove B molecule
  I4 cell(-2);
  Mol.Set( mix->pM, Cidx,cell ); 

  // Set velocity of third body (Molecule M)
  Reac->InnerProb->vX=Reac->InnerProb->Vm3 - Vr*P2->c_mass[P1->id];
  Mol.Set(Reac->InnerProb->pX , mix->Vidx, Reac->InnerProb->vX );
}//end cEIntDiscr::RealizeRecombination

