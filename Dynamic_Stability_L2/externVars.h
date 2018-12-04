#ifndef __externVars_h__
#define __externVars_h__

#include "DataTypes.h"
#include "defs.h"

extern struct SSpecies stSpecies[];
extern S32Bit	n_C, n_N;
extern FLOAT8	C[Large], N[Large], B[Large], PhiC[Large], PhiN[Large];
extern S32Bit	nSpecies,Counter, TotalVars;
extern FLOAT8	dMinBiomassThreshold, dProbCon, IntStr;
extern FLOAT8	dDiffConst;
extern FLOAT8	dFlux1, dFlux2;
extern FLOAT8	R_0;
extern FLOAT8	Affinity11, Affinity12, Affinity21, Affinity22;
extern FLOAT8	Y11, Y12, Y21, Y22, Alpha1, Alpha2;
extern S32Bit	idum, TotalVars;
extern S32Bit	nMaxIter;
extern S32Bit	nPNvalue;
extern FLOAT8 	TotalTime;
extern S32Bit	nChemicalSetSize;	
extern S32Bit	nSecreteSetSize;
extern S32Bit	nAffectingSetSize;
extern FLOAT8	fChemMutationProb;	//Mutation Prob.of SecreteSet & AffectingSet
extern FLOAT4	ran1();
extern FILE	*fSpeciesData, *fResourceData, *fSpeciesMatrix, *fNetwork, *fAfterPerturb, *fTotalNoOfStSt, *fFinalStats;
#endif // __externVars_h__
