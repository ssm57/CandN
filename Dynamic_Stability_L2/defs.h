#ifndef __defs_h__
#define __defs_h__

#include "DataTypes.h"

#define	Large						1000			// max num of species

//#define n_C							5
//#define n_N							5


#define dBIOMASS_MAX			0.0009
#define dBIOMASS_MIN			0.0001

#define	C_MAX			1.0
#define	C_MIN			1.0

#define	N_MAX			1.0
#define	N_MIN			1.0

#define	B_MAX			1e-3
#define	B_MIN			1e-3

#define	dLambda_MAX			11.0
#define	dLambda_MIN			10.0

#define	MonodN_MAX				0.1
#define	MonodN_MIN				0.01

#define	MonodC_MAX				0.1
#define	MonodC_MIN				0.01

#define	dMIN_BIOMASS_THRESHOLD	5.0e-3

#define dMINYDOTTHRESHOLD		1.0e-12

#define dMAXTIME				2.0e10

#define dDELTA_MAX	 			0.0009
#define dDELTA_MIN	 			0.0004

#define fNEWBIOMASSTIMES		10.0

#define	dOneIterationTime		0.001
/*
#define delta_birth	 			0.06
#define delta_death  			0.06

*/
//------- for CVODE
# define RTOL 	RCONST(1.0e-11) 	// scalar relative tolerance
# define ATOL 	RCONST(1.0e-15) 	// vector absolute tolerance components
# define T0 	RCONST(0.0) 		// initial time
# define T1 	RCONST(100.0) 		// first output time
# define TMULT RCONST(10.0) 		// output time factor
# define NOUT 	5 						// number of output time

/*
typedef struct SSpecies
{
	FLOAT8	dBirthRate;
	FLOAT8	dDeathRate;
	FLOAT8	dBiomass;
	FLOAT8	dChem;
} ss;
*/
typedef struct SSpecies
{
	//FLOAT8	dYield[N+1];//Yield[0] is the Yield on the R, Yield[1] is Yield on C1
	FLOAT8	LambdaC;
	FLOAT8	LambdaN;
	FLOAT8	YC;	
	FLOAT8	YN;	
	FLOAT8	dBiomass;
	char	Presence;
} ss;

#endif // __defs_h__
