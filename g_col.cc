/*
 g_col.cc
 Collision class for DSMC C++ project 
 25.05.98	Kashkovsky A.V.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

#include "mytype.h"
//#include "mathd.h"
//#include "tmatrix.h"
#include "reader.h"
//#include "dump.h"
//#include "array.h"

#include "g_gas.h"
#include "particles.h"
#include "cells.h"
#include "g_col.h"
#include "g_mesh.h"
#include "conf.h"
#include "g_index.h"
#include "level2.h"

#include "g_move.h"
#include "flowf.h"
#include "rndm.h"
#include "Paralleler/paral_parent.h"

extern cgGases    *Gases; // Gases list
extern cConf      *Conf;  // Configuration of the calculation
extern cParticles  Mol;   // Particles
extern cCells      Cells; // Cells
extern cgMesh     *Mesh;  // Calculation mesh
extern cgIndex    *Index; // Indexation
extern cgLev2     *Lev2;  //  Subcell (only for collision)
extern cRndm Rn;
extern cgMove *Move; // Move
extern cgFlowField *Field;  // Flowfields mesh
extern ParalParent *paralleler;     // Paralleler father


//********************************************** cgCol

/********************************************************
 * cgCol::cgCol						*
 *  Constructor						*
 * 10.05.00			Kashkovsky A.V.		*
 ********************************************************/
cgCol::cgCol()
{
  read = new cRead();
  SetCollVel( "VHS" );
  read->AddStr( "MolModel", CollVelName, 31);
  mm = 0;
// gas has to be read before!!!!
  ncomp = Gases->Num();
  ptc1L=ptc1M=ptc2L=ptc2M=NULL;
  int nn = (1+ncomp)*ncomp/2; // lenght of thriangle matrix
  iCrp = Cells.AddMember("CRP",'D',8*nn,"CRP ");
  Vidx = Mol.AddMemberV8( "V", "Vel " );
  len_subcell=0;
  subcell=0;
  N2gas=-1;
}

/********************************************************
 * cgCol::SetCollVel					*
 *  Set Collision veloscity calculation function						*
 * 10.05.00			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::SetCollVel(
 const char *funcName	// name of the function
)
{
  if( funcName == NULL || funcName[0] == 0 ||
          !strcmp(funcName, "VHS" ) )
  {colvel=&cgCol::CollVelVHS;
   strcpy( CollVelName, "VHS" );
   mm = 0;
  }
  else if(!strcmp(funcName, "VSS" ) )
  {colvel=&cgCol::CollVelVSS;
   strcpy( CollVelName, "VSS" );
   mm = 1;
  }
  else
  {ostringstream oss;
   oss << "Unknow Collision velocity calculation function: "
       << funcName;
   OnErr(1,"cgCol::SetCollVel",oss.str().c_str() );
  }
}

/********************************************************
 * cgCol::Make						*
 *  Prepere all collision parameters			*
 * 10.05.00			Kashkovsky A.V.		*
 * 01.02.03			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::Make( R8 Temp)
{
  SetCollVel( CollVelName );
  R8 val = pow( 100., 1-2*Gases->Get(0)->alpha);
  for( int i=0;i<Cells.Num();++i)
  {  SetCrp( i, &val ); 
  }
}

/********************************************************
 * cgCol::SetCrp
 * cgCol::GetCrp
 *  set/get max realative velocity Crp for ONE component
 * 16.06.10   Kashkovsky A.V.
 ********************************************************/
void  cgCol::SetCrp( const int &cell, R8 *val)
{ Cells.Set( cell, iCrp, *val );
}
void  cgCol::GetCrp( const int &cell, R8 *val)
{ Cells.Get( cell, iCrp, val );
}

/********************************************************
 * cgCol::Collide					*
 *  General Collision of the particles without internal *
 *  energy and chemestry				*
 * 06.04.00			Kashkovsky A.V.		*
 * 06.06.10			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::Collide(
  R8 Step		// time step
)
{

  I4 lev(0);       // level of 2-nd mesh
  R8 ncc;
  R8 np;        // number of particles in cell
  R8 np2;       // np2=np-1 - for second ptc gener.
  R8 Cr;        // Relative powered velocity
  R8 Crp;       // Relative powered velocity in the cell
  R8 ncol;      // number of coll. chekings
  R8 tau;       // sub step for collision selection
  R8 vol;       // volume of the cell

  gL=gM=Gases->Get(0);		// gas
  nCol = 0;

// Background cell loop 
//  for( bgcell=Cells.Num()-1;bgcell>=0;--bgcell)
  for( bgcell=0;bgcell<Cells.Num();++bgcell)
  {
    const I4 glob_cell = paralleler->global_number_from_local( bgcell, paralleler->GetIProc());
    subcells = Lev2->GetCellLevel( bgcell );
    if( subcells > 1 ) // bg cell split
    { // reindex
      lev =  Lev2->GetCellLevel( bgcell, 0 );
      SubIndex( bgcell );
    }
    for(cell=0; cell<subcells;++cell)
    {
      if( subcells > 1 )
      { vol = Mesh->SubVol( glob_cell, cell, 1, lev);
        np = SubNumPtc( cell );
      }
      else
      { vol = Mesh->Vol( glob_cell );
        np = Index->NumPtc( bgcell );
      }
      if(np < 2. ) continue;
      // constant for calculation coll. checking number
      ncc= Conf->FNum* gL->c_STr[0] * np*(np-1)/2. /
           vol * Step;
// for Weight scheme
      R8 rweight   = Mesh->WeightRad(glob_cell);
      ncc*=rweight;
// after adding of paralleler
      np-=1.e-5; np2=np-1.;


      // the calculation number of collisions
/*************************************
 Number of collision from majorant frequency
$$
 Mf = N*(N-1)/2 * Fn * (\sigma_T*Cr)_{max}*\Delta t/ Vol
$$

 where                                    \\
  $N$  - number of model particle in cell  \\
  $Fn$ - ratio real particles to model     \\
  $\Delta t$ - time step                   \\
  $Vol$ - Volume of the cell               \\
$\sigma_T$ - realative cross section       \\
  $Cr$ - relativ velocity                  \\


 from: Bird (equvation (4.60))
$$
 \sigma_T = \pi d^2
$$
Bird (equvation (4.63))
$$
 d = d_{ref} \sqrt{ (2 k T_{ref}/(mr*Cr^2))^\alpha / 
                     \Gamma(2-\alpha ) }
$$
 extract references parameter, which don't changes in 
calculation:
$$
 STr = \pi d_{ref}^2 *  (2 k T_{ref}/(mr))^\alpha / 
                     \Gamma(2-\alpha )
$$
$\sigma_T$ will be:
$$
 \sigma_T    = STr*Cr^{(-2*\alpha)}
$$
or
$$
 \sigma_T*Cr = STr*Cr^{(1-2*\alpha)}
$$
or, by replacemen $Crp = Cr^{(1-2*\alpha)}$

$$
 \sigma_T*Cr = STr*Crp
$$

$ (\sigma_T*Cr)_{max} $ found and store for each background cell
  
*************************/
      GetCrp( bgcell, &Crp );
      ncol = ncc*Crp;
int real_coll, fict_coll;
real_coll=fict_coll=0;
      // loop for checking collisions
      for(tau=-log(Rn.dm()); tau<ncol; tau-=log(Rn.dm()))
      { fict_coll++;
      // search two candidates for collision

        oL = (I4)(np*Rn.dm());
        oM = (I4)(np2*Rn.dm());
        if( oM>=oL ) oM++;
        PretendPtc(oL,0,pL,ptc1L,ptc2L);
        PretendPtc(oM,0,pM,ptc1M,ptc2M);
     // calculate relative and mean velocity
        Mol.Get( pL, Vidx, &vL );
        Mol.Get( pM, Vidx, &vM );
        Vrel = vL - vM;
        Vmod2 = Vrel*Vrel;
        if( Vmod2 < 1.e-15 ) continue;
        Vmod  = sqrt(Vmod2);
        Vmean = vL*gL->c_mass[gM->id]
              + vM*gM->c_mass[gL->id];
        Cr = pow( Vmod, 1-2*gL->alpha );
        if( Cr > Crp )
        { Crp = Cr;
          SetCrp( bgcell, &Crp );
          ncol = ncc*Crp;
        }
        else if( Cr/Crp < Rn.dm() ) continue; // no coll.

        // the are collision
        nCol++;
real_coll++;
       if (Field) Field->SampColl(bgcell,gL->id,gM->id);

        // Calculate new velocity of the particles
        CollVel();
      } // for( tau=-log(Rn.dm()); tau<mf; tau-=log(Rn.dm()) )
/*
if( bgcell > 140 && bgcell < 150 )
{ printf("g_coll:cell %d nptc %g crp %g ncol %g fict %d real %d \n",
bgcell, np, Crp,
ncol, fict_coll, real_coll );
}
*/

    } // end for(cell=0; cell<subcells;++cell)
  } // end for( bgcell=Cells.Num()-1;bgcell>=0;--bgcell)
  Conf->CountStep();  // for protocol prn from cColl, Chem, ions..
}

/********************************************************
 * cgCol::CollVelVHS						*
 *  Variable Hard Spheres collision			*
 * 06.05.00			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::CollVelVHS()
{
  R8 A,B,C;
  V8 vCr;

//  B is the cosine of a random elevation angle
  B=2.*Rn.dm()-1.;
  A=sqrt(1.-B*B);
  vCr.x  = B*Vmod;
//  C is a random azimuth angle
  C=2.*M_PI*Rn.dm();
  vCr.y = A*cos(C)*Vmod;
  vCr.z = A*sin(C)*Vmod;
  vL = Vmean + vCr*gM->c_mass[gL->id];
  vM = Vmean - vCr*gL->c_mass[gM->id];
  Mol.Set( ptc2L, Vidx, vL );
  Mol.Set( ptc2M, Vidx, vM );
}

/********************************************************
 * cgCol::CollVelVSS					*
 *  Variable Soft Spheres collision			*
 * 06.05.00			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::CollVelVSS()
{
  R8 A,B,C,D,cs,sn;
  V8 vCr;

//--B is the cosine of the deflection angle for the VSS model, eqn (11.8)
  B=2.*pow( Rn.dm(), 1./gL->c_vss[gM->id] )-1.;
  A=sqrt(1.-B*B);
  C=2.*M_PI*Rn.dm();
  cs=cos(C);
  sn=sin(C);
  D=sqrt( Vrel.y*Vrel.y + Vrel.z*Vrel.z);
  if( D > 1.e-6 )
  {
    vCr.x=B*Vrel.x+A*sn*D;
    vCr.y=B*Vrel.y+A*(Vmod*Vrel.z*cs-Vrel.x*Vrel.y*sn)/D;
    vCr.z=B*Vrel.z-A*(Vmod*Vrel.y*cs+Vrel.x*Vrel.z*sn)/D;
  }
  else
  {
    vCr.x = B*Vmod;
    vCr.y = A*cs*Vmod;
    vCr.z = A*sn*Vmod;
  }
  vL = Vmean + vCr*gM->c_mass[gL->id];
  vM = Vmean - vCr*gL->c_mass[gM->id];
#warning  TO DO: write to pointer in VSS
  Mol.Set( pM, Vidx, vM );
  Mol.Set( pL, Vidx, vL );
}

// * * * * * * * * * * * * * * * * Second level indexation


/********************************************************
 * cgCol::SubIndex					*
 *  Flash ptc address					*
 * 06.02.01			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::SubIndex(
 I4 bgcell	// BG cell for make 2-nd level index
)
{
  int subcells;  // total number of cells in this bg
  subcells = Lev2->GetCellLevel( bgcell );
  if( subcells > len_subcell )
  { len_subcell = subcells;
    if( subcell ) delete [] subcell;
    try{ subcell = new AI4[ len_subcell*ncomp ]; }
    catch(std::bad_alloc)
    {ostringstream oss;
     oss << "Can't allocate " 
         << len_subcell*ncomp
         << " bytes memory";
     OnErr(1,"cgCol::SubIndex",oss.str().c_str() );
    }
  }
  else 
  {
    I4 i,nsub=subcells*ncomp;
    for(i=0;i<nsub;i++)subcell[i].Num()=0;
  }
// get particles from BG index 
  Mesh->InitFor2LevIndex( bgcell );
  I4 i,c,pid,sub;
  I4 nmol;
  V8 pos;
  for( c=0;c<ncomp;c++)
  { nmol=Index->NumPtc(bgcell,c);
    for( i=0;i<nmol;++i)
    { pid=Index->PtcAdr( bgcell, i, c );
      pos = Move->PtcCoor( pid );
      sub = Mesh->Cell2( &pos, bgcell);
if( sub < 0 )
{
printf("proc: %d g_col.cc:369, bgcell = %d, sub = %d < 0, subcells %d\n",
 paralleler->GetIProc(), bgcell, sub, subcells);

int realcell, proccell, global_cell,wh_is;
I4 Cidx  = Mol.IndexOf( "Cell");
Mol.Get( pid, Cidx, &realcell );
    proccell = paralleler->
          local_number_from_global(realcell);
  global_cell = paralleler->global_number_from_local(
      proccell,   paralleler->GetIProc() );
  wh_is=paralleler->whose_is_cell(global_cell);

printf("mol cell %d, local cell %d global_cell %d proc %d whose_is_cell %d\n",
 realcell, proccell, global_cell, paralleler->GetIProc(), wh_is );

 Mol.Print( pid );
  exit(2);
}


      SubAdd( sub, pid, c );
    }
  }
}

/********************************************************
 * cgCol::SubAdd						*
 *  Add particle adress for given component					*
 * 06.02.01			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::SubAdd(
 I4 cell,	// sub cell number
 I4 ptc,	// particle ID
 I4 nc	// component number
)
{ I4 rc, i = cell*ncomp + nc;
  rc = subcell[i].Add(ptc);
  if( rc == SPACE )
  {ostringstream oss;
   oss << "Error in add new adress cell " 
       << ios_base::dec << cell
       << "for component "
       << ios_base::dec << nc;
   OnErr(1,"cgCol::SubAdd",oss.str().c_str() );
  }
}

/********************************************************
 * cgCol::SubDel						*
 *  Del particle adress for given component					*
 * 04.06.02			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::SubDel(
 I4 cell,	// sub cell number
 I4 ptc,	// particle address
 I4 nc	// component number
)
{ I4 i = cell*ncomp + nc;
  subcell[i].Del(ptc);
}

/********************************************************
 * cgCol::SubNumPtc						*
 *  Get number of particle for given component	*
 * 08.02.01			Kashkovsky A.V.		*
 ********************************************************/
I4 cgCol::SubNumPtc(
 I4 cell,	// subcell 
 I4 nc	// component number
)
{ I4 i =  cell*ncomp + nc;
  return(subcell[i].Num());
}

/********************************************************
 * cgCol::SubPtc						*
 *  Get particle number for given component and count	*
 * 07.02.01			Kashkovsky A.V.		*
 ********************************************************/
I4 cgCol::SubPtc(
 I4 cell,	// subcell
 I4 count,	// count number of particle 
 I4 nc	// component number
)
  { I4 i = cell*ncomp + nc;
  return(subcell[i].Get(count));
}



// * * * * * * * * * * * * * * * * Dump/Restore functions


//********************************************** cgCol

/********************************************************
 * cgCol::Dump
 *  Dump collisions parameters
 * 29.01.02			Kashkovsky A.V.		*
 ********************************************************/
void cgCol::Dump(
 const char *file,	// file name for dump
 I4 Step, // number of full steps made for dump time
 R8 Time  // Modeling time for dump time
)
{
/***-------------------- TODO
  const char *Func="cgCol::Dump";
  FILE *fp;
  int ver=0, i;

  OPEN_DUMP( file )
  OUT_HEAD("COLL", ver )

  OUT_NUM( "Step", Step )
  OUT_NUM( "Time", Time )
  OUT_NUM( "Ncel", ndcl )

  OUT_STRUCTI( ndcl )

  OUT_EOF

-------------------- TODO ________***/
}


/********************************************************
 * cgCol::Restore
 *  Restore collision parameters
 * 29.01.02			Kashkovsky A.V.
 ********************************************************/
void cgCol::Restore(
 const char *file,	// file for write
 I4 *Step,	// current step
 R8 *Time	// current time
)
{
/***-------------------- TODO
  const char *Func="cgCol::Restore";
  FILE *fp;
  int rc;
  const char key[5]={0};
  I4 Ver,len,nnn=-1;

  OPEN_REST(file, "COLL", Ver)

  int i;

// loop for reading all keywords
// Some date dump only for information, and we don't need 
// to use them and we make duumy read
  while(1)
  {
    RED_KEY( key, len )
    if( !strncmp( key, "ENDF", 4 ) ) break;
    else GET_NUM( "Step", Step)
    else GET_NUM( "Time", Time)
    else GET_NUM( "Ncel", &nnn )
    else GET_STRUCTI( nnn, ndcl ) // try to restore array members
  }
  fclose(fp);
-------------------- TODO ________***/
}



/********************************************************
 * cgCol::PretendPtc
 *  define variables of one ptc of the collisional pair
 * maintains of the pL , ptc1L ,  ptc2L variables
 * igL & igM are already defined here
 * 29.01.17			 Shevyrin Al.An.
 ********************************************************/
void cgCol::PretendPtc(const I4 &o,     // number in cell array
		       const I4 &gi,    // gas index for molecule
		       I4 &p_set,       // set to index in Mol
		       PP &ptc1_set,    // ptc pointer for read
		       PP &ptc2_set    // ptc pointer for write
		       )
{
  //  cout<<"------>Pretend o "<< o <<" gi "<<gi;
  p_set = ( subcells > 1 ? SubPtc( cell, o, gi ) 
	    : Index->PtcAdr( bgcell, o, gi ) );
  ptc2_set = ptc1_set = Mol.pMolGet(p_set);
//   cout<<" <------- Pretend p "<< p_set <<" ptc "
//       <<(void*)ptc_set<<endl;
}// end cgCol::PretendPtc

/********************************************************
 * cgCol::ReversePair
 * swipe L&M molecules in dissociation
 * 29.01.17			 Shevyrin Al.An.
 ********************************************************/
void cgCol::ReversePair(){
  I4 dummy;
  dummy = oL;  oL = oM;   oM = dummy;
  dummy = pL;  pL = pM;   pM = dummy;
  dummy = igL; igL = igM; igM = dummy;
  pcgGas pg;
  pg = gL;  gL = gM; gM = pg;
  PP ptc;
  ptc = ptc1L ; ptc1L = ptc1M ; ptc1M = ptc;
  ptc = ptc2L ; ptc2L = ptc2M ; ptc2M = ptc;
}

/*************************************************
 unity vector of random direction 25.1.17. Shevyrin Al.An.
*************************************************/
void UnitVecXYZ(V8 &r){
  R8 C,S,Theta;
  C = 1.-2.*Rn.dm();
  S = sqrt(1.-C*C);
  Theta = 2.*M_PI*Rn.dm();
  r.x = C;
  r.y = S*cos(Theta);
  r.z = S*sin(Theta);
}; // end UnitVecXYZ
