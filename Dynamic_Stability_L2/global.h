#ifndef __global_h__
#define __global_h__

#include "defs.h"

struct	SSpecies stSpecies[Large] = {0};
S16Bit	n_C, n_N;
FLOAT8	C[Large], N[Large], B[Large], PhiC[Large], PhiN[Large];
S32Bit	nSpecies, Counter, TotalVars;
FLOAT8	dFlux1, dFlux2;
FLOAT8	dMinBiomassThreshold;
FLOAT8  TotalTime;
FLOAT8	dProbCon, R_0, dDiffConst,ResourceFlux;
FLOAT8	Affinity11, Affinity12, Affinity21, Affinity22;
FLOAT8	Y11, Y12, Y21, Y22, Alpha1, Alpha2;
S32Bit	nSecreteSetSize, TotalVars;
S32Bit	nAffectingSetSize;
FLOAT8	fChemMutationProb;
S32Bit	idum;
S32Bit	nMaxIter;
S32Bit	nPNvalue;
FILE	*fSpeciesData, *fResourceData, *fSpeciesMatrix, *fNetwork, *fAfterPerturb, *fTotalNoOfStSt, *fFinalStats;
extern float ran1();

#endif // __global_h__
