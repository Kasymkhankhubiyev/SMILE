/*
 g_col.h
 general collision classes for DSMC C++ project 
 10.05.00           Kashkovsky A.V.
 07.06.10 version 2 Kashkovsky A.V.
*/


#ifndef G_COL_H
#define G_COL_H

typedef char*  PP;

/*************************************************
 general collisions 
*************************************************/
class cgCol
{ public:
  cRead *read;	// reader
  cStruct dump;

  I4 mm;	// molecular model (0-VHS, 1-VSS)
  I4 iCrp;	// Crp index in cell container
  I4 Vidx;      // index of ptc velosity in container
  I4 nCol;	// Number of collisions per step

// collisions variable for prticle L & M
  I4 oL;	// Order number of particle in cell list
  I4 oM;
  I4 pL;	// Index of ptc. in particle container
  I4 pM;
  I4 igL,igM;
  pcgGas gL;	// Gas description
  pcgGas gM;

  //Shevyrin 26.12.2016 for improved version
//collision pair (before collision molecules)
  PP ptc1L;  //ptr to oL 
  PP ptc1M;  //ptr to oM
//collision pair (after collision molecules)
  PP ptc2L;  //ptr to L - molecule after coll 
  PP ptc2M;  //ptr to M - molecule after coll 

  // save N2 gas index for ionization Tv(N2), -1 if no N2
  I4 N2gas; 

  //Shevyrin 27.1.2017 
  // define variables of one ptc of the collisional pair
  virtual void PretendPtc(const I4 &o,     // number in cell array
			  const I4 &gi,    // gas index for this molecule 
			  I4 &p_set,       // set to index in Mol
			  PP &ptc1_set,    // ptc pointer for read
			  PP &ptc2_set);   // ptc pointer for write
  virtual void ReversePair(); // swipe L&M molecules in dissociation

  V8 vL;       // velocity of particle
  V8 vM;

// general collision parameters
  I4 bgcell;    // current BG cell
  I4 subcells;  // number of subcells
  I4 cell;      // current 2-nd level cell
  R8 vol;       // volume of the cell


  V8 Vmean;	// Mean mass velocity
  V8 Vrel;	// Relative velocity
  R8 Vmod;	// Module  of Relative velocity
  R8 Vmod2;	// Module$^2$ of Relative velocity

// second level indexation
  I4 ncomp;     // number of gass component
  AI4 *subcell; // particle numbers in sub cell
  I4   len_subcell; // length of subcell

  cgCol();
  virtual ~cgCol(){delete read; delete [] subcell;};
  virtual void Make(R8 Temp);


// collision velocity functions
  char CollVelName[32];	// name of the collision functions
  void (cgCol::*colvel)();

  virtual void SetCollVel(const char *Name );
  virtual void CollVel()
       {return( (this->*colvel)() );};

  virtual void CollVelVHS();
  virtual void CollVelVSS();

// Max relative velosity for componets Crp
  virtual void SetCrp( const int &cell, R8 *val);
  virtual void GetCrp( const int &cell, R8 *val);


// 2-nd level indexation
  virtual void SubIndex(I4 BGcell);
  virtual void SubAdd(I4 sub, I4 ptc, I4 ncomp=0);
  virtual void SubDel(I4 sub, I4 ord, I4 ncomp=0);
  virtual I4 SubNumPtc(I4 sub, I4 c=0);
  virtual I4 SubPtc(I4 sub, I4 count=0, I4 c=0 );

// collision
  virtual void UpdateInt( I4 i){};
  virtual void Collide(R8 Step);

  virtual void Dump(const char *file,I4 Step,R8 Time);
  virtual void Restore(const char *file,I4 *Step,R8 *Time);


};

/*************************************************
 unity vector of random direction 25.1.17. Shevyrin Al.An.
*************************************************/
void UnitVecXYZ(V8 &r);

#endif
