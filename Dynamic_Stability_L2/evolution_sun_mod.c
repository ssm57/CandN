//-------- includes for CVODE
#include <cvode/cvode.h>   
#include <nvector/nvector_serial.h>  // serial N_Vector types, fcts., macros 
#include <cvode/cvode_dense.h>       // prototype for CVDense 
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM 
#include <sundials/sundials_types.h> // definition of type realtype

#define Ith(v, i)    NV_Ith_S(v, i-1)     // Ith numbers components 1..nSpecies 
#define IJth(A, i, j) DENSE_ELEM(A, i-1, j-1) // IJth numbers rows, cols 1..nSpecies


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "externVars.h"

// ------------- for random number generator
#define IA 		16807
#define IM 		2147483647
#define AM 		(1.0/IM)
#define IQ 		127773
#define IR 		2836
#define NTAB 	32
#define NDIV 	(1+(IM-1)/NTAB)
#define EPS 	1.2e-7
#define RNMX 	(1.0-EPS)

/* some global variables */
realtype dMaxydot;

//------- for CVODE
realtype	reltol;
realtype	t;
static	N_Vector	y=NULL;
static	N_Vector	ydot=NULL;
static	N_Vector	abstol=NULL;
static	void		*cvode_mem=NULL;
int	flag;
S32Bit nCheckMinydot;
FLOAT8	Cdot[Large], Ndot[Large], Bdot[Large];


static int check_flag(void *flagvalue, char *funcname, int opt);

/* ---------------------------------------------------- */
/* random number generator, returns a float between 0 and 1 */
/* for random init idum should be initialised to a negative integer */

//-------------------------------------------------------------------------------------------------------------------//
//											ran1
//-------------------------------------------------------------------------------------------------------------------//

float ran1(long *idum)
{
	S32Bit 	j;
	long 	k;
	float temp;
	static long iy=0;
	static long iv[NTAB];

	if (*idum<=0 || !iy)
	{
		if ( -(*idum) < 1 ) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7; j>=0; j--)
		{
			k = (*idum) / IQ;
			*idum = IA * (*idum-k * IQ) - IR * k;
			if (*idum<0) *idum += IM;
			if (j<NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}

	k = (*idum)/IQ;
	*idum = IA * (*idum - k*IQ) - IR*k;
	if (*idum<0) *idum += IM;
	j = iy/NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ( (temp = AM*iy) > RNMX ) return RNMX;
	else return temp;
}

double GaussRand(FLOAT8 mu, FLOAT8 sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      //printf("\n half value would have been %e but %e chosen\n", mu, mu + sigma * (double) X1);
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;
  return (mu + sigma * (double) X1);
}

//-------------------------------------------------------------------------------------------------------------------//										
//											nPopulationRates
//-------------------------------------------------------------------------------------------------------------------//

static int nPopulationRates(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	int i, c, n, b;
	FLOAT8 Min;
	dMaxydot=0.0;
	// -------------------------------- making all the negetive biomass equal to 0.0

	for (i=1; i<= TotalVars; i++)	
	{	if ( Ith(y, i) < 1e-20 )
		{	Ith(y, i) = 0.0;
			continue;
		}
	}

	//Ith(ydot, 0) =	0.0;

	//Copy all the Ith y
	for(c=1; c<=n_C; c++)
		C[c] = Ith(y,c);
	for(n=1; n<=n_N; n++)
		N[n] = Ith(y, n_C + n);
	for(b=1; b<= n_N*n_C; b++)
		B[b] = Ith(y, n_C + n_N + b);

	for(c=1; c<= n_C; c++)
		Cdot[c] = PhiC[c] - (dDiffConst*C[c]);
	for(n=1; n<= n_N; n++)
		Ndot[n] = PhiN[n] - (dDiffConst*N[n]);
	
	b = 0;
	for(c=1; c<= n_C; c++)
	{	for(n=1; n<= n_N; n++)
		{	b++;
			if ( (stSpecies[b].LambdaC * C[c]) < (stSpecies[b].LambdaN * N[n]) )
				Min = stSpecies[b].LambdaC * C[c];
			else
				Min = stSpecies[b].LambdaN * N[n];
				
			Bdot[b] = B[b] * ( Min - dDiffConst );
			Cdot[c] -= ( (double)(1.0/stSpecies[b].YC) * B[b] * Min );
			Ndot[n] -= ( (double)(1.0/stSpecies[b].YN) * B[b] * Min );			
		}
	}

	for(c=1; c<=n_C; c++)
		Ith(ydot,c) = Cdot[c];
	for(n=1; n<=n_N; n++)
		Ith(ydot, n_C + n) = Ndot[n];
	for(b=1; b<= n_N*n_C; b++)
		Ith(ydot, n_C + n_N + b) = Bdot[b];
	
	if(nCheckMinydot == 1)
	{	dMaxydot=fabs(Ith(ydot, 1));
		for(i=2; i<= TotalVars; i++)
		{
			if(fabs(Ith(ydot, i))>dMaxydot)
			dMaxydot=fabs(Ith(ydot, i));		
		}

		if(dMaxydot <= dMINYDOTTHRESHOLD)
		nCheckMinydot = 0;
		//printf("\n**** ydot has gone below threshold at t=%f****", t);
	}

	return (0);
}

//-------------------------------------------------------------------------------------------------------------------//
//						vDyneig
//-------------------------------------------------------------------------------------------------------------------//

void vDyneig() 
{
	S32Bit	i, c, n, b;
	//FLOAT8 Temp1, Temp2;
	static	realtype th;
	U32Bit nMaxSteps = 1e8, nSteps;

	y = N_VNew_Serial( TotalVars);
	//TotalVars is n_N + n_C + (n_N*n_C);

	if (check_flag((void *)y, "N_VNew_Serial", 0)) exit(1);

	abstol = N_VNew_Serial(TotalVars);

	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) exit(1);

	// Initialize y
	//Ith(y, 0) = 0.0;

	for(c=1; c<=n_C; c++)
		Ith(y,c) = C[c] ;
	for(n=1; n<=n_N; n++)
		Ith(y, n_C + n) = N[n];
	for(b=1; b<= n_N*n_C; b++)
		Ith(y, n_C + n_N + b) = stSpecies[b].dBiomass;

	nCheckMinydot = 1;	// TO check when to stop integrating by checking maximum ydot.

	// Set the scalar relative tolerance 
	reltol = RTOL;

	// Set the vector absolute tolerance
	for(i=1; i <= TotalVars; i++)
		Ith(abstol , i) = ATOL ;

	// Create solver memory and specify BDF and Newton Iteration 
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) exit(1);

	// Initialize the integrator memory & specify derivative function
	flag = CVodeInit(cvode_mem, nPopulationRates, T0, y);

	if (check_flag(&flag, "CVodeInit", 1)) exit(1);

	// Specify scalar relative tolerance & vector absolute tolerances
	flag = CVodeSVtolerances(cvode_mem, reltol, abstol);

	if (check_flag(&flag, "CVodeSVtolerances", 1)) exit(1);

	// Specify the CVDENSE dense linear solver
	flag = CVDense(cvode_mem, TotalVars);

	if (check_flag(&flag, "CVDense", 1)) exit(1);

	// Specify maxsteps
	flag = CVodeSetMaxNumSteps(cvode_mem, 2e9);

	// Specify maxnef
	flag = CVodeSetMaxErrTestFails(cvode_mem, 1e9);

	// Specify maxncf
	flag = CVodeSetMaxConvFails(cvode_mem, 1e9);

	//------------------------------------------------------------------------------------------------------------------------
	// ---------------------------------- integrating ------------------------------------------------------------------------
	t=0.0;
	th=0.0; nSteps = 0;

	/*
	if(counter == 3)
	{	printf("	t			C1			C2			N1			N2			B1			B2			B3			B4\n");
		printf("%.5e	%.5e	%.5e	%.5e	%.5e	%.5e	%.5e	%.5e	%.5e\n", 
		t, Ith(y, 1), Ith(y, 2), Ith(y, n_C+1), Ith(y, n_C+2), Ith(y, n_C+n_N+1), Ith(y, n_C+n_N+2), Ith(y, n_C+n_N+3), Ith(y, n_C+n_N+4));
	}*/	
	
	while(nCheckMinydot == 1 && nSteps <= nMaxSteps)
	{	//SEE PAGE 28 OF cv.guide2.pdf about the 2nd and 4th arguments.
		flag = CVode(cvode_mem, th + dOneIterationTime, y, &t, CV_ONE_STEP);

		if (check_flag(&flag, "CVode", 1)) 
		{
			fprintf( stderr, "error integrating\n" );
			exit(1);
		}
		
		/*
		if(counter == 3)	
		{	printf("%.5e	%.5e	%.5e	%.5e	%.5e	%.5e	%.5e	%.5e	%.5e\n", 
			t, Ith(y, 1), Ith(y, 2), Ith(y, n_C+1), Ith(y, n_C+2), Ith(y, n_C+n_N+1), Ith(y, n_C+n_N+2), Ith(y, n_C+n_N+3), Ith(y, n_C+n_N+4));}*/
		
		if(t > dMAXTIME)
		{	nCheckMinydot = 0;
			//printf("\nydots deceasing very slowly & Time of integration has exceeded the dMAXTIME,"
			//" so we stop at this instant of time\n");
			//printf(" dMaxydot = %e\n", dMaxydot);
		}

		th += dOneIterationTime;
		nSteps += 1;
	}

	TotalTime = t;
	if (nSteps >= nMaxSteps)
		printf("\nOOPS\n");

	for(c=1; c<=n_C; c++)
		C[c] = Ith(y,c);
	for(n=1; n<=n_N; n++)
		N[n] = Ith(y, n_C + n);
	for(b=1; b<= n_N*n_C; b++)
		stSpecies[b].dBiomass = Ith(y, n_C + n_N + b);

	//printf("\ndone in t=%e, nSteps=%ld\n", t, nSteps);
	
	//printf("	y[0]= %.15f\n\n", Ith(y, 0));
	//for (i=1; i<= nSpecies; i++) 
		//printf("	y[%ld]= %.15f\n", i, Ith(y, i));


	/* Free y and abstol vectors */
	N_VDestroy_Serial(y);
	N_VDestroy_Serial(abstol);

	/* Free integrator memory */
	CVodeFree(&cvode_mem);

}

static int check_flag(void *flagvalue, char *funcname, int opt)
{
	int *errflag;

	// Check if SUNDIALS function returned NULL pointer - no memory allocated 
	if (opt == 0 && flagvalue == NULL) 
	{
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
		return(1); 
	}
	else // Check if flag < 0 
	{
		if (opt == 1) 
		{
			errflag = (int *) flagvalue;
			if (*errflag < 0) 
			{
				fprintf(stderr, "\n SUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
				return(1); 
			}
		}
		else 
		{
			if (opt == 2 && flagvalue == NULL) 
			// Check if function returned NULL pointer - no memory allocated 
			{
				fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
				return(1); 
			}
		}
	}
	return(0);
}
