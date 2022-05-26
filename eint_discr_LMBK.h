/* eint_discr_LMBK.h
 *  VV12 exchange, SSH model
 *  28.06.2016 Molchanova A.N.
 *  17.03.2020 Latyshev A.K.
 */

#ifndef EINT_DISCR_LMBK_H
#define EINT_DISCR_LMBK_H


#include "eint_discr.h"

class cEIntDiscrLMBK: public cEIntDiscr
{ public:

 cEIntDiscrLMBK();


 virtual void EintExch(cgCol *col); 	// returns 1 if there was exchange, 0 if not
 void ColVV12( // make VV12-exchange
			 pcgGas gA,	// gas property for A particle - VV12 exchange
			 pcgGas gB,	// gas property for B particle - catalyst
			 U1 *Lv,	// Levels of vibrational energy [number of vib modes] of A particle
			 R8 &Erel,	// Relative translational energy
			 int increment // -1/+1 => decrease/increase level of first mode: "-1" => i1-1 , i2+2; "+1" => i1+1 , i2-2
			);
 void Split_lv2(		// divide total level of bending mode on two levels of submodes uniformly
			   int lv2,		 // total level of bending mode (lv2 = lv2_1+lv2_2)
			   int maxlevel, // maximal number of level (in first or in second - the same)
			   U1 &lv2_1,	 // level of first submode (define)
			   U1 &lv2_2	 // level of second submode (define)
			  );
 R8 Prob_VV12_decrease(	// calculate probability of VV12: lv1->lv1-1, lv2->lv2+2
			   pcgGas gA,	 // gas property for A particle (VV12-exchange)
			   pcgGas gB,	 // gas property for B particle (catalyst)
			   R8 Erel,      // relative translational energy
			   int lv1,		 // energy level of symmetric mode (0 mode) - decrease by 1
			   int lv2		 // total energy level of bending mode (sum of 1st and 2nd) - increase by 2
			  );
 R8 Prob_VV12_increase(	// calculate probability of VV12: lv1->lv1+1, lv2->lv2-2
			   pcgGas gA,	 // gas property for A particle (VV12-exchange)
			   pcgGas gB,	 // gas property for B particle (catalyst)
			   R8 Erel,      // relative translational energy
			   int lv1,		 // energy level of symmetric mode (0 mode) - increase by 1
			   int lv2		 // total energy level of bending mode (sum of 1st and 2nd) - decrease by 2
			  );

};
#endif
