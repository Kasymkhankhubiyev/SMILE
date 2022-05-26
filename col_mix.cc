/*
 mix.cc
 Collision class for DSMC C++ project 
 20.03.02	Kashkovsky A.V.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

#include "mytype.h"
//#include "mathd.h"
#include "reader.h"
//#include "dump.h"
//#include "array.h"
//#include "tmatrix.h"
#include "g_gas.h"
#include "particles.h"
#include "cells.h"
//#include "g_col.h"
#include "col_mix.h"
#include "g_eint.h"
#include "g_mesh.h"
#include "conf.h"
#include "gasweight.h"
#include "g_chem.h"
#include "g_index.h"
#include "level2.h"
#include "flowf.h"
#include "rndm.h"
#include "Paralleler/paral_parent.h"

extern cgGases *Gases;	// Gases list
extern cgReacs *Reacs;	// Reaction list
extern cConf   *Conf;	// Configuration of the calculation
extern cParticles  Mol;   // Particles
extern cCells      Cells; // Cells
extern cgMesh     *Mesh;  // Calculation mesh
extern cgIndex    *Index; // Indexation
extern cgLev2     *Lev2;  //  Subcell (only for collision)
extern cgEInt *EIntMod; // model of internal energy
extern cGasWeight *gasWeights; // Weights for trace components
extern cReacsOfDirIonization *ReacsDirIon; // Direct ionization
extern cReacsOfIonRecomb *Reacs_iRecomb; // Dissociative recombination
extern cgMeanTv    *meanTv; // Counter mean Tv in cells
extern cgFlowField *Field;  // Flowfields mesh
extern cTe_Tlocal *Te_local;  // 



extern cRndm Rn;

#ifdef TAG
// Debug tracing of special particle life 
extern I8 debug_tag; // tag for printing particle history
#endif


extern ParalParent *paralleler;     // Paralleler father
//********************************************** cMix

/********************************************************
 * cMix::Make						*
 *  Prepere all collision parameters			*
 * 20.03.02			Kashkovsky A.V.		*
 ********************************************************/
void cMix::Make(R8 Temp)
{ 
/*
  cgCol::Make( Temp);
  !!! we can't call this function, because we don't want to
  create 'crp' array for single species
  So, we need to copy content of cgCol::Make.
*/

  int i;
  SetCollVel( CollVelName );

// create other crp storage
/* in cgCol::cgCol was create Crp store for
  (ncomp+1)*ncomp/2 array. we need to fill it
  Create single array with a value and copy for each cell
*/
  int nn = (ncomp+1)*ncomp/2;
  R8 *acrp = new R8[ nn ];
  R8 val = pow( 100., 1-2*Gases->Get(0)->alpha);
  for(i=0;i<nn;++i) acrp[i]=val;
  for(i=0;i<Cells.Num();++i)
  {  Cells.SetStr( i, iCrp, (char *)acrp ); 
  }
  delete [] acrp;

  mf = new R8[ncomp*ncomp];
  npm = new R8[ncomp];
  nColM=new I4[ncomp*ncomp];
  for( i=0;i<ncomp;i++)nColM[i] = 0;
  ChemFactor.Init(ncomp);
  ChemFactor.Set(1.);

  N2gas=0;
  while( strcmp(Gases->Get(N2gas)-> name ,"N2") ) {
    ++N2gas;
    if (N2gas>=Gases->Num()) {
      N2gas=-1; // no N2 gas
      break;
    }
  }
}

/********************************************************
 * cMix::SetCrp
 * cMix::GetCrp
 *  set/get max realative velocity Crp for ONE component
 * 16.06.10   Kashkovsky A.V.
 ********************************************************/
void  cMix::SetCrp( const int &cell, R8 *val)
{ R8 *p;
  p = (R8*)Cells.PGet( cell, iCrp );
// index in the half of thriangle matrix
  int i, row,col;
  if( gL->id <= gM->id ) {col = gL->id; row = gM->id;}
  else                   {row = gL->id; col = gM->id;}
  i = row*(row+1)/2+col;

  p[i] = *val;
}
void  cMix::GetCrp( const int &cell, R8 *val)
{  R8 *p;
  p = (R8*)Cells.PGet( cell, iCrp );
// index in the half of thriangle matrix
  int i, row,col;
  if( gL->id <= gM->id ) {col = gL->id; row = gM->id;}
  else                   {row = gL->id; col = gM->id;}
  i = row*(row+1)/2+col;

  (*val) = p[i];
}


/********************************************************
 * cMix::Collide					*
 *  Collision of the mixture particles with internal
 *  energy and chemestry				*
 * 20.03.03			Kashkovsky A.V.		*
 * 16.06.10			Kashkovsky A.V.		*
 ********************************************************/
void cMix::Collide(
  R8 Step		// time step
)
{ 

  R8 Cr;        // Realtive powered velocity
  R8 Crp;       // Realtive powered velocity in the cell
  R8 ncol;      // number of coll. chekings
  R8 tau;       // sub step for collision selection
  R8 ff;
  int i,j(0);
  I4 lev(0);       // level of 2-nd mesh
  int change(0);

  nCol = 0;

  if( Reacs ) {Reacs->TestingReset();  Reacs->nChemReac=0;}
  if(gasWeights) gasWeights->BucketReset();
//  for( i=0;i<ncomp*ncomp;i++)nColM[i] = 0;

// Background cell loop 
//  for( bgcell=Cells.Num()-1;bgcell>=0;--bgcell)
  for( bgcell=0;bgcell < Cells.Num();++bgcell)
  {// update probability (for T-depended ZrZv )
//     printf("*****bcell %5d **",bgcell);
      {
          // Check, is there at least 1 particle in the cell or not
          bool particleIs(false);
          for( int iComp=0; iComp < ncomp; ++iComp) {
              if ( Index->NumPtc( bgcell, iComp ) > 0 ) {
                  particleIs = true;
                  break;
              }
          }
          if ( ! particleIs ) continue;
      }
    const I4 glob_cell = paralleler->global_number_from_local( bgcell, paralleler->GetIProc());

    if(Gases->VarZrZv) EIntMod->ProbUpdate(bgcell);

    subcells = Lev2->GetCellLevel( bgcell );
    if( subcells > 1 ) // bg cell split
    { // reindex
      lev =  Lev2->GetCellLevel( bgcell, 0 );
      SubIndex( bgcell );
    }
    for(cell=0; cell<subcells;++cell)
    {
      if( subcells > 1 )
	vol = Mesh->SubVol(glob_cell, cell, 1, lev );
      else
	vol = Mesh->Vol( glob_cell);
      FixCellPtcNum();
      if(np < 2. ) continue;
      Majorant( vol, Step, ncol);

// for Weight scheme
      R8 rweight   = Mesh->WeightRad(glob_cell);
      ncol*=rweight;
//
// debug
int fict_coll, real_coll;
fict_coll=real_coll=0;
    // loop for checking collisions
      for(tau=-log(Rn.dm()); tau<ncol; tau-=log(Rn.dm()))
      {
fict_coll++;
      // search a component for collision
        ff = Rn.dm();
        for( i=0; i<ncomp; i++ )
        { for( j=0; j<=i; j++ )
          { if( ff <= mf[i+j*ncomp] ) break;
          }
          if( j<=i ) break;
        }
	igL=i; igM=j;
      // search two candidates for collision
        if( i==j )
        { if( npm[i] < 2. ) continue;
          oL = (I4)( (npm[i]-1.e-6)   *Rn.dm());
          oM = (I4)( (npm[j]-1.-1.e-6)*Rn.dm());
          if( oM>=oL ) oM++;
        }
        else
        {
          oL = (I4)( (npm[i]-1.e-6)   *Rn.dm());
          oM = (I4)( (npm[j]-1.e-6)   *Rn.dm());
        }
	PretendPtc(oL,i,pL,ptc1L,ptc2L);
	gL=Gases->Get(i);
	PretendPtc(oM,j,pM,ptc1M,ptc2M);
	gM=Gases->Get(j);

#ifdef TAG
       Mol.PrintTag(pL, debug_tag, "Collide, pair pL" ,
              paralleler->GetIProc() );
       Mol.PrintTag(pM, debug_tag, "Collide, pair pL" ,
              paralleler->GetIProc() );
#endif

       // sample Ev for TCE
       if ( meanTv) {
	 R8 Er,Ev;
	 EIntMod->GetE(pL,&Er,&Ev ); meanTv->CountEv(bgcell,i,Ev);
	 EIntMod->GetE(pM,&Er,&Ev ); meanTv->CountEv(bgcell,j,Ev);
       }

      // calculate relative and mean velocity
        Mol.Get( ptc1L, Vidx, &vL );
        Mol.Get( ptc1M, Vidx, &vM );
	Vrel = vL - vM;
        Vmod2 = Vrel*Vrel;
        // sample T, if it need  
        EIntMod->AddMeanTTrn(bgcell, this ); 

        if( Vmod2 < 1.e-15 ) continue;
        Vmod  = sqrt(Vmod2);
        Vmean = vL*gL->c_mass[gM->id]
              + vM*gM->c_mass[gL->id];

        Cr = pow( Vmod, 1-2*gL->c_alpha[gM->id] );
        GetCrp( bgcell, &Crp ); // !!! Crp for gL+gM component
        if( Cr > Crp )
        { Crp = Cr;
          SetCrp( bgcell, &Crp );
          change=1; // !!! need to recalculate Majorant !!!
        }
        else if( Cr/Crp < Rn.dm() ) continue;  // no coll.
        // the are collision
// debug
real_coll++;
        nCol++;
        nColM[gL->id+gM->id*ncomp]++;
       if (Field) Field->SampColl(bgcell,gL->id,gM->id);

      // check chemical reaction

        if( Reacs )
        {
#ifdef TAG
       Mol.PrintTag(pL, debug_tag, "Collide, before chem pL" ,
              paralleler->GetIProc() );
       Mol.PrintTag(pM, debug_tag, "Collide, before chem pM" ,
              paralleler->GetIProc() );
#endif

          if( Reacs->Chemistry( this ) )
          { // there was a reaction,
            // recalculate number of
            // particles by component
            // because it was changed
	    FixCellPtcNum();
            //rebuild the majorant fr.
            // and continue collision 
            Majorant( vol, Step, ncol);
            ncol*=rweight;	// for weight scheme
            continue;
          }
        } // end  if( Reacs )
#ifdef TAG
       Mol.PrintTag(pL, debug_tag, "Collide, after chem pL" ,
              paralleler->GetIProc() );
       Mol.PrintTag(pM, debug_tag, "Collide, after chem pM" ,
              paralleler->GetIProc() );
#endif
      // set r/w pointers in gasweights for collision
       if (gasWeights) gasWeights->Train(this);
	
      // make exchange of Internal energy

	if( Rn.dm()*ChemFactor.Data(gL->id,gM->id)>1. ) {
	  --nCol;
          --nColM[gL->id+gM->id*ncomp];
	  continue;
	}
	
	EIntMod->EintExch( this );
#ifdef TAG
       Mol.PrintTag(pL, debug_tag, "Collide, EintExch pL" ,
              paralleler->GetIProc() );
       Mol.PrintTag(pM, debug_tag, "Collide, EintExch pM" ,
              paralleler->GetIProc() );
#endif

      // Calculate new velocity of the particles
	CollVel();
        if( change )
        {  Majorant( vol, Step, ncol);
           ncol*=rweight;	// for weight scheme
           change=0;
        }
#ifdef TAG
       Mol.PrintTag(pL, debug_tag, "Collide, CollVel pL" ,
              paralleler->GetIProc() );
       Mol.PrintTag(pM, debug_tag, "Collide, CollVel pM" ,
              paralleler->GetIProc() );
#endif

      } // for( tau=-log(Rn.dm()); tau<mf; tau-=log(Rn.dm()) )

// dbug
/*
if( glob_cell > 140 && glob_cell < 150 )
{ printf("col_mix,cell %d nptc %g crp %g ncol %g fict %d real %d \n",
glob_cell, np, Crp,
ncol, fict_coll, real_coll );
}
*/

   if (gasWeights)  gasWeights->Push2Mol();
   if (Te_local) Te_local->SetTe_to( meanTv->TvCell(bgcell,N2gas) ); //Tv_N2
   if(ReacsDirIon) ReacsDirIon->Ionization(this,Step,Conf->FNum/vol*rweight);
   if(Reacs_iRecomb) Reacs_iRecomb->IonsDissRecombination(this,Step,Conf->FNum/vol*rweight);

    } // end for(cell=0; cell<subcells;++cell)
  } // End for(mesh->curcell=0; mesh->curcell<ncell; mesh->curcell++)
  if( meanTv ) meanTv->ResetTvCell();

//if( Reacs ) Reacs->TestingPrint();
  Conf->CountStep();  // for protocol prn from cColl, Chem, ions..
}

/********************************************************
 * cMix::FixCellPtcNum
 * charge npm[] array, calc np
 * 17.1.17 Shevyrin AlAn
 ********************************************************/
void cMix::FixCellPtcNum(){
  I4 i;
  if( subcells > 1 )
    { for( np=i=0;i<ncomp;++i)
	{ npm[i]=SubNumPtc( cell,i );
	  np+=npm[i];
	}
    }
  else
    { for( np=i=0;i<ncomp;++i)
	{ npm[i]=Index->NumPtc( bgcell,i );
	  np+=npm[i];
	}
    }
}  

/********************************************************
 * cMix::FixCellPtcNum
 * charge npm[comp_] array
 * np become invalid
 * 22.10.18 Shevyrin AlAn
 ********************************************************/
void cMix::UpdCellPtcNum_comp(const I4 &comp_){
  if( subcells > 1 )
    {
      npm[comp_]=SubNumPtc( cell,comp_);
    }
  else
    {
      npm[comp_]=Index->NumPtc( bgcell,comp_ );
    }
}  


/********************************************************
 * cMix::Majorant
 *  Calculation of the majorant fr. and probabilityes
 * 20.03.03			Kashkovsky A.V.		*
 ********************************************************/
void cMix::Majorant(
 const R8 &vol,    // volume of the cell
 const R8 &Step,   // time step
       R8 &ncol    // 1/Mf
)
{
  I4 i,j;
  R8 Cr;		// Realtive powered velocity
  R8 ff;		// Realtive powered velocity


/***********************
 for each c pairs found max. possible number of collisions
 Number of collision from majorant frequency
 if( i == j ) 
$$ 
 Mf_{ij} = N_i*(N_j-1)/2 * Fn * (\sigma_{Tij}*Cr)_{max}*\Delta t/ Vol
$$
 if( i != j )
$$ 
 Mf_{ij} = N_i*N_j * Fn * (\sigma_{Tij}*Cr)_{max}*\Delta t/ Vol
$$
(look g\_col.cc: cgCol::Collide for addition explanation
***********************/
 ncol = 0.;
 for( i=0; i<ncomp; i++ )
 {gL = Gases->Get(i);
  if( npm[i] < 1. ) 
  { for( j=0; j<=i; j++ ) mf[j+i*ncomp]=mf[i+j*ncomp]=0.;
    continue;
  }
  for( j=0; j<=i; j++ )
  { gM = Gases->Get(j);
    GetCrp( bgcell, &Cr ); // !!! Crp for gL+gM component
                           // max.rel.vel^1-2.\alpha
//     Cr=500.;

//     if ((i==j)&&(npm[i]==100.0)) {
//       printf("PPPP npm %d\n Fn %e\n all %e\n chf %e STr %e Cr %e \n",

// (I4)npm[i],Conf->FNum,
// 	     ChemFactor.Data(i,j)* gL->c_STr[i] /
// 	     vol * Step * Cr,
// 	     ChemFactor.Data(i,j),gL->c_STr[i],Cr);
//     }

    if( i == j )
    { 
      mf[i+i*ncomp]= ChemFactor.Data(i,j)*Conf->FNum* gL->c_STr[i] * npm[i]*(npm[i]-1.)/2. /
                    vol * Step * Cr;
    }
    else 
    { 
      mf[j+i*ncomp]=
      mf[i+j*ncomp]= ChemFactor.Data(i,j)*Conf->FNum
	* gL->c_STr[j] * npm[i]*npm[j] 
	/ vol * Step * Cr;
    }
    if (gasWeights){
      mf[i+j*ncomp]=mf[j+i*ncomp]*=gasWeights->CollKoef(j,i);
    }
    ncol+=mf[i+j*ncomp];
  }
 }
 if( ncol < 1.e-9 ) {ncol=0.;return;}

// make a histogram from number of collisions
 ff=0.;
 for( i=0; i<ncomp; i++ )
 { for( j=0; j<=i; j++ )
   {
     ff+=mf[i+j*ncomp]/ncol;
     mf[i+j*ncomp] = ff;
   }
 }

 i=j=ncomp-1; 
 mf[i+j*ncomp]=1.;
}

/********************************************************
 * cMix::DelPtcInCurCell
 *  Delete particle from index list in current cell 
 * 17.06.10        Kashkovsky A.V.
 ********************************************************/
void  cMix::DelPtcInCurCell( 
 const I4 &ptc, // order index of particle in the cell
 const I4 &gas  // index of gas component
)
{
  if( subcells > 1 ) // 2-nd level
  { SubDel( cell, ptc, gas );
  }
  else  // 1-st level
  { Index->Del(bgcell, ptc, gas);
  }
  npm[gas]--;
}

/********************************************************
 * cMix::AddPtcInCurCell
 *  Add particle index to index list in current cell 
 * 17.06.10        Kashkovsky A.V.
 ********************************************************/
void  cMix::AddPtcInCurCell( 
 const I4 &ptc, // global index of particle in the cell
 const I4 &gas  // index of gas component
)
{
  if( subcells > 1 ) // 2-nd level
  { SubAdd( cell, ptc, gas );
  }
  else  // 1-st level
  { Index->Add(bgcell, ptc, gas);
  }
  npm[gas]++;
}
