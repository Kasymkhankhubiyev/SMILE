/*
 col_mix.h (mix.h)
 mixture collision classes for DSMC C++ project 
 20.3.02   Kashkovsky A.V.
 09.06.10   Kashkovsky A.V.
*/


#ifndef COL_MIX_H
#define COL_MIX_H

#include "g_col.h"
#include "tmatrix.h"

/******************************************************
Mixture collisional 
******************************************************/
class cMix: public cgCol
{ public:
//  I4 nmix == ncomp;      // number of mixtures
  R8 *mf;       // number of coll. for component  
  I4 *nColM;	// Number of collisions per step for component

  //Shevyrin 26.12.2016 for improved version 

  //NB:   w/o gas weights
  //ptc2L=ptc1L and ptc2M=ptc1M
  //write ptc2LM after final read!
  //with weights ptc2-s are buffers for resulting mols
  //

  // np,npm are defined and valid in scope of cMix::Collide
  R8 *npm;      // number of particles for component
  R8 np;        // number of particles in cell
// becouse of chemistry, it is need to del/add
// particles in curent cell. It is done by 
// DelPtcInCurCell and AddPtcInCurCell
// but it is need to know level of cell


  cMix(){nColM=NULL;mf=npm=NULL;np=0.;};
  virtual ~cMix(){ delete [] nColM;
                   delete [] mf; 
                   delete [] npm;
                 };

// collision
  virtual void Collide(R8 Step);

  void FixCellPtcNum(); // charge npm[] array,calc np
  void UpdCellPtcNum_comp(const I4 &comp_); // charge npm[comp_] array, not np 

  virtual void Majorant(
    const R8 &vol, 
    const R8 &Step,
          R8 &ncol);

// for chemistry. change gas components in cell
  virtual void DelPtcInCurCell( const I4 &ptc, const I4 &gas ); 
  virtual void AddPtcInCurCell( const I4 &ptc, const I4 &gas ); 

// Max relative velosity for componets Crp
  virtual void SetCrp( const int &cell,  R8 *val);
  virtual void GetCrp( const int &cell,  R8 *val);

  virtual void Make( R8 Temp);

// may be, it will be neeed
//  virtual void AddTemp(I4 cell,R8 T);


  
// 14.03.11                Shevyrin Al.An.
// next functions togeather with variable ChemFactor
// inroduced to deal with probr > 1
//   wich appear when tabulated crossection used.
//   basic idea: 
//  1. if cgReacTabDiss defined the matrix ChemFactor={1} does exist.
//  2. if reaction probr > 1 then 
//  A. ChemFactor*= for current components redifened by factor probr,
//  B. c_STr[i]=c_STr[j]*=probr (dummy increase in crossection), 
//     which means inrease of ncol and mf[] in function Collide here
//  C. probr in Prob() authomatically decreased
//  D. Before realization of real collision in Collide()
//     make rejection of 1/ChemFactor collisions
  tTmatrix<R8> ChemFactor;
};
#endif

