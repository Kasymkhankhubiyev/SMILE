/* eint_discr_LMBK.cc
 *  VV12 exchange, SSH model
 *  28.06.2016 Molchanova A.N.
 *  17.03.2020 Latyshev A.K
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mytype.h"
#include "mathd.h"
#include "reader.h"
//#include "array.h"
//#include "tmatrix.h"
#include "g_gas.h"
#include "g_chem.h"
#include "col_mix.h"
#include "cells.h"
#include "particles.h"
#include "rndm.h"
#include "eint_discr_LMBK.h"

extern cgGases *Gases;      // Gases list
extern cParticles   Mol;    // Particles
extern cRndm Rn;     		// random number generator
extern R8 Boltz; // Boltzmann constant (k)

using namespace std;
// data to calculate VV12 probabilities
const double A_1(1./2.);
const double A_2(8./11.);
const double Masredvib_1(1.46e-26);	// reduced mass of symmetric mode
const double Masredvib_2(1.338e-26);	// reduced mass of bending mode
const double Nu_1(3.855e+13);			// frequency of symmetric mode
const double Nu_2(2.0e+13);				// frequency of bending mode
const double r0(3.763e-10);				// Lennard-Jones potential parameter, need to calculate next alpha
const double alpha(17.5/r0);			// repulsive potential parameter
const double planck(6.626e-34);         // Planck's constant

/*******************************************************
 *  cEIntDiscrLMBK::cEIntDiscrLMBK
 *  constructor
 *  17.03.20       Latyshev A.K.
 *******************************************************/
cEIntDiscrLMBK::cEIntDiscrLMBK():cEIntDiscr()
{
    Counter_VV12 = 0;
//  Counter_VV23_increase_total = 0;
//  Counter_VV23_decrease_0 = 0;
//  Counter_VV23_increase_0 = 0;
//  Counter_VV23_decrease_1 = 0;
//  Counter_VV23_increase_1 = 0;
//  Counter_VV23_decrease_2 = 0;
//  Counter_VV23_increase_2 = 0;
//  Counter_VV23_decrease_3 = 0;
//  Counter_VV23_increase_3 = 0;
//  Counter_VV23_decrease_4 = 0;
//  Counter_VV23_increase_4 = 0;
//  Counter_VV23_decrease_5 = 0;
//  Counter_VV23_increase_5 = 0;
}

/*******************************************************
 * cEIntDiscrLMBK::EintExch
 * Exchange of internal energy
 * 17.03.20       Latyshev A. K.
 *******************************************************/
void cEIntDiscrLMBK::EintExch(
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

  U2 Lr_L[2], Lr_M[2];	// Rotational energy levels for L and M particles
  U1 Lv_L[8], Lv_M[8];	// Vibrational energy levels for L and M particles

  Mol.GetStr( col->ptc1L, LVidx, (char *)Lv_L );
  Mol.GetStr( col->ptc1M, LVidx, (char *)Lv_M );

  Erel = col->gL->c_massr[col->igM]*col->Vmod2*0.5;


  // Calculate standart RT probability
  R8 Prob_RT_L = col->gL->rotmods*ProbRT.Data( col->gL->id, col->gM->id);
  R8 Prob_RT_M = col->gM->rotmods*ProbRT.Data( col->gM->id, col->gL->id);

  // Calculate VV12 probability
  R8 Prob_VV12_L[2];
  R8 Prob_VV12_M[2];

  Prob_VV12_L[0] = Prob_VV12_decrease(col->gL,col->gM,Erel,Lv_L[0],Lv_L[1]);
  Prob_VV12_L[1] = Prob_VV12_increase(col->gL,col->gM,Erel,Lv_L[0],Lv_L[1]);
  Prob_VV12_M[0] = Prob_VV12_decrease(col->gM,col->gL,Erel,Lv_M[0],Lv_M[1]);
  Prob_VV12_M[1] = Prob_VV12_increase(col->gM,col->gL,Erel,Lv_M[0],Lv_M[1]);

  if (Prob_VV12_L[0]>1) Counter_VV12++;
  if (Prob_VV12_L[1]>1) Counter_VV12++;
  if (Prob_VV12_M[0]>1) Counter_VV12++;
  if (Prob_VV12_M[1]>1) Counter_VV12++;
  
  int increment[2]={-1,1}; // to define decrease/increase level of mode: in VV12

  // Calculate total probability of internal energy exchange
  R8 P_EintExch_total(0);
  // Add probabilities of RT- and VT-exchange
  P_EintExch_total+=ProbSum.Data(col->igL, col->igM);
  // Add probabilities of VV12-exchange
  P_EintExch_total+=Prob_VV12_L[0]+Prob_VV12_L[1];
  P_EintExch_total+=Prob_VV12_M[0]+Prob_VV12_M[1];

 if( col->gL->type == ATOM && col->gM->type == ATOM ) return; 	// there wasn't internal energy exchange - go to usual elastic collision

// check if there is internal energy transfer according to LB scheme
 if ( Rn.dm() > P_EintExch_total ) return; 	// there wasn't internal energy exchange - go to usual elastic collision

// check with first molecule
 Rnd = Rn.dm();
 Dust = Prob_RT_L/P_EintExch_total;

 // RT transfer for L molecule
 if( Rnd < Dust )
 {
   Mol.GetStr( col->ptc1L, LRidx, (char *)Lr );
   ColRot(col->gL,col->gM, Lr, Erel );
   Mol.SetStr( col->ptc2L, LRidx, (char *)Lr );
   goto END;
 }

 // VT transfer for L molecule
 Rnd-=Dust;
 for( m=0;m<col->gL->vibmods;m++)
 { Dust = ProbVT[m].Data(col->gL->id, col->gM->id)/
          P_EintExch_total;
   if( Rnd < Dust )
   {
     Mol.GetStr( col->ptc1L, LVidx, (char *)Lv );
     ColVib(col->gL,col->gM, m,Lv, Erel );
     Mol.SetStr( col->ptc2L, LVidx, (char *)Lv );
     goto END;
   }
   Rnd-=Dust;
 } // end for( m=0;m<col->gL->vibmods;m++)

  // Make VV12-exchange
  // VV12 exchange for L molecule
  for( int i=0;i<2;i++)
  { // decrease/increase level cycle
    Dust = Prob_VV12_L[i] / P_EintExch_total;
    if( Rnd < Dust )
     { // exchange
		Mol.GetStr( col->ptc1L, LVidx, (char *)Lv);
        ColVV12(col->gL, col->gM, Lv, Erel, increment[i] );
        Mol.SetStr( col->ptc2L, LVidx, (char *)Lv);
        //Counter_VV12++;
        goto END;
      }
    Rnd-=Dust;
  }

// check with second molecule
 Dust = Prob_RT_M/P_EintExch_total;

 // RT transfer for M molecule
 if( Rnd < Dust )
 {
   Mol.GetStr( col->ptc1M, LRidx, (char *)Lr );
   ColRot(col->gM,col->gL, Lr, Erel );
   Mol.SetStr( col->ptc2M, LRidx, (char *)Lr );
   goto END;
 }

 // VT transfer for M molecule
 Rnd-=Dust;
 for( m=0;m<col->gM->vibmods;m++)
 { Dust = ProbVT[m].Data(col->igM, col->igL)/
          P_EintExch_total;
   if( Rnd < Dust )
   {
     Mol.GetStr( col->ptc1M, LVidx, (char *)Lv );
     ColVib(col->gM,col->gL, m,Lv, Erel );
     Mol.SetStr( col->ptc2M, LVidx, (char *)Lv );
     goto END;
   }
   Rnd-=Dust;
 } // end for( m=0;m<col->gM->vibmods;m++)

  // VV12 exchange for M molecule
  for( int i=0;i<2;i++)
  { // decrease/increase level cycle
   Dust = Prob_VV12_M[i] / P_EintExch_total;
   if( Rnd < Dust )
     { // exchange
	   Mol.GetStr( col->ptc1M, LVidx, (char *)Lv );
       ColVV12(col->gM, col->gL, Lv, Erel, increment[i]);
       Mol.SetStr( col->ptc2M, LVidx, (char *)Lv );
       //Counter_VV12++;
       goto END;
     }
   Rnd-=Dust;
  }

END:
// change relative velocity

 W = sqrt(2.* Erel/ col->gL->c_massr[col->igM]);
 F = W/col->Vmod;
 col->Vrel*=F;
 col->Vmod =W;
}


/*******************************************************
 *  cEIntDiscrLMBK::ColVV12
 *  Make VV12 exchange
 *  17.03.2020			Latyshev A. K.
 *******************************************************/
void cEIntDiscrLMBK::ColVV12(
							pcgGas gA,		// gas property for A particle - VV12 exchange
							pcgGas gB,		// gas property for B particle - catalyst
							U1 *Lv,			// Levels of vibrational energy [number of vib modes] - for A particle
							R8 &Erel,		// Relative energy
							int increment	// flag -1/+1 => decrease/increase level of first mode: "-1" => i1-1 , i2+2; "+1" => i1+1 , i2-2
						   )
{

  // total energy for split over Trans. and Vib
  R8 Ecol = Erel + gA->EV[0][ Lv[0] ] + gA->EV[1][ Lv[1] ];

  Lv[0]+=increment;    // -1/+1
  Lv[1]-=increment*2;    // +2/-2

  // divide new level on two modes uniformly
  //Split_lv2(lv2,gA->LevelV[1]-1,Lv[1],Lv[2]);

  //recount relative translational energy
  Erel = Ecol - (gA->EV[0][ Lv[0] ] + gA->EV[1][ Lv[1] ]);
}


/*******************************************************
 *  cEIntDiscrLMBK::Split_lv2
 *  Divide total level of bending mode on two levels of submodes uniformly
 *  29.06.2016			Molchanova A.N.
 *******************************************************/
void cEIntDiscrLMBK::Split_lv2(		// divide total level of bending mode on two levels of submodes uniformly
							  int lv2,		 // total level of bending mode (lv2 = lv2_1+lv2_2)
							  int maxlevel,	 // maximal number of level in submode (in first or in second - the same)
							  U1 &lv2_1,	 // level of first submode (need to define)
							  U1 &lv2_2	     // level of second submode (need to define)
							 )
{
  if( lv2 <= maxlevel )
    {
	  lv2_1 = (int)(((double)(lv2+1))*Rn.dm());
    }
  else  // define region of possible values of first submode
    {
      R8 diff = lv2-maxlevel;
      R8 reg = maxlevel - diff + 1;
	  lv2_1 = (int)(reg*Rn.dm() + diff);
    }
  lv2_2 = lv2 - lv2_1 ;
}




/*******************************************************
 *  cEIntDiscrLMBK::Prob_VV12_decrease
 *  Calculate probability of VV12, decrease: lv1->lv1-1, lv2->lv2+2
 *  17.03.2020			Latyshev A.K.
 *******************************************************/
R8 cEIntDiscrLMBK::Prob_VV12_decrease(	// calculate probability of VV12: lv1->lv1-1, lv2->lv2+2
									 pcgGas gA,	 	// gas property for A particle (VV12-exchange)
									 pcgGas gB,	 	// gas property for B particle (catalyst)
									 R8 Erel,     	// relative translational energy
									 int lv1,		// energy level of symmetric mode (0 mode) - decrease by 1
									 int lv2		// total energy level of bending mode (sum of 1st anf 2nd) - increase by 2
									)
{
  R8 prob(0.);
  int lv2_1,lv2_2;
//  if (strcmp(gA->name,"CO2")==0) return 0;
  if ( (lv1 < 1) || (lv2 > (gA->LevelV[1]-3))) return 0;
  if (Erel< (2*(gA->EV[1][1]-gA->EV[1][0]) - (gA->EV[0][1]-gA->EV[0][0]))) return 0;
  //prob = pow(alpha*A_2/(M_PI),4)*pow(A_1/(Masredvib_2*Nu_2),2)*planck*gA->c_massr[gB->id]*Erel/(256*Masredvib_1*Nu_1);
  prob = pow(alpha/(M_PI),6)*pow(A_2,4)*pow(A_1/(Masredvib_2*Nu_2),2)*pow(planck,3)/(2048.*Masredvib_1*Nu_1);
  lv2_1 = (int)(((double)(lv2+1))*Rn.dm());
  lv2_2 = lv2 - lv2_1;
  prob = prob * lv1*((lv2_1+2)*(lv2_1+1)+4*(lv2_1+1)*(lv2_2+1)+(lv2_2+2)*(lv2_2+1));
  //prob = prob * lv1*(lv2+1)*(lv2+2)*(lv2+3); 
  //prob = prob * lv1*(lv2+1)*(lv2+2)/300.; 
  return prob;
}


/*******************************************************
 *  cEIntDiscrLMBK::Prob_VV12_increase
 *  Calculate probability of VV12, increase: lv1->lv1+1, lv2->lv2-2
 *  17.03.2020			Latyshev A.K.
 *******************************************************/
R8 cEIntDiscrLMBK::Prob_VV12_increase(	// calculate probability of VV12: lv1->lv1+1, lv2->lv2-2
									 pcgGas gA,		// gas property for A particle (VV12-exchange)
									 pcgGas gB,		// gas property for B particle (catalyst)
									 R8 Erel,      	// relative translational energy
									 int lv1,		// energy level of symmetric mode (0 mode) - increase by 1
									 int lv2		// total energy level of bending mode (sum of 1st anf 2nd) - decrease by 2
									)
{
  R8 prob(0.);
  int lv2_1,lv2_2;
 // if (strcmp(gA->name,"CO2")==0) return 0;
  if ( (lv2 < 2) || (lv1 > (gA->LevelV[0]-2))) return 0;
  if (Erel< ( (gA->EV[0][1]-gA->EV[0][0]) - 2*(gA->EV[1][1]-gA->EV[1][0]) ) ) return 0;
  //prob = pow(alpha*A_2/(M_PI),4)*pow(A_1/(Masredvib_2*Nu_2),2)*planck*gA->c_massr[gB->id]*Erel/(256*Masredvib_1*Nu_1);
  //prob = prob * (lv1+1)*lv2*(lv2-1)/300.;
  prob = pow(alpha/(M_PI),6)*pow(A_2,4)*pow(A_1/(Masredvib_2*Nu_2),2)*pow(planck,3)/(2048.*Masredvib_1*Nu_1);
  lv2_1 = (int)(((double)(lv2+1))*Rn.dm());
  lv2_2 = lv2 - lv2_1;
  prob = prob*(lv1+1)*(lv2_1*(lv2_1-1)+4*lv2_1*lv2_2+lv2_2*(lv2_2-1));
  return prob;
}
