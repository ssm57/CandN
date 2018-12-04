#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "global.h"

extern void vDyneig();
extern void vKillSpecies(S32Bit, S32Bit*);
extern void vAddSpecies(S32Bit, S32Bit*, char*);
extern void	Perturbation();

extern S32Bit nCheckNewChem(S32Bit, S32Bit, S32Bit, S32Bit*, S32Bit**, S32Bit, S32Bit, S32Bit*);
extern double GaussRand (FLOAT8, FLOAT8);
char U[Large][3]; char V[Large][3];
S16Bit G[Large][3] = {0};//G[b][1] is the Carbon that species-b uses and G[b][2] is the Nitrogen it uses
S32Bit	NoOfSwitching[10];
	
/* --------------------- MAIN FUNCTION ------------------------ */

int main(int argc, char *argv[])
{
	S32Bit	c, n, b, NoOfPlus=0, State;


	idum	=				-atol( argv[1]);//Random Seed
	n_C		=				atoi(argv[2]);//No of Carbon sources
	n_N		=				atoi(argv[3]);//No of Nitrogen sources

	dDiffConst = 			1.0;//atof(argv[4]);//The Diffusion constant

	//printf("\nMaxNoOfExp=%ld	TotalNoOfRandOrderExp=%ld",	MaxNoOfExp,	TotalNoOfRandOrderExp);
	//printf("\nSteady state is reached when no. of plus becomes Zero\n");
	
	TotalVars = n_N + n_C + (n_N*n_C);

	char	s6[5000] = "FinalStats";
	sprintf(s6, "%s_%s_%s_%s.txt", s6, argv[1], argv[2], argv[3]);
	fFinalStats = fopen(s6, "w");

	FLOAT8	Min = 0.0;

	b=0;
	for (c=1; c<=n_C; c++)
		for (n=1; n<=n_C; n++)
		{	b++;
			G[b][1] = c, G[b][2] = n;
		}

			//***********************Set the parameters	**************************************
			//Set the Lambdas and Y

			stSpecies[1].LambdaC=41.0,	stSpecies[1].LambdaN=16.0, stSpecies[1].YC=	0.37, stSpecies[1].YN=	0.27;
			stSpecies[2].LambdaC=35.0,	stSpecies[2].LambdaN=50.0, stSpecies[2].YC=	0.64, stSpecies[2].YN=	0.10;
			stSpecies[3].LambdaC=52.0,	stSpecies[3].LambdaN=27.0, stSpecies[3].YC=	0.47, stSpecies[3].YN=	0.22;
			stSpecies[4].LambdaC=56.0,	stSpecies[4].LambdaN=44.0, stSpecies[4].YC=	0.14, stSpecies[4].YN=	0.59;			


			//PhiC[1]=263.1905368268, PhiC[2]=762.4358329177, PhiN[1]=753.6025048494, PhiN[2]=909.2989107966, State=33;
			PhiC[1]=42.379017, PhiC[2]=767.297924, PhiN[1]=939.785674, PhiN[2]=154.161797, State=32;
			//PhiC[1]=725.686599, PhiC[2]=513.760442, PhiN[1]=890.846736, PhiN[2]=748.544332, State=31;
			//PhiC[1]=1.594447, PhiC[2]=389.434418, PhiN[1]=696.546566, PhiN[2]=12.304911, State=30;
			//PhiC[1]=494.482706, PhiC[2]=51.033900, PhiN[1]=464.981384, PhiN[2]=527.401854, State=29;
			//PhiC[1]=930.506055, PhiC[2]=384.118564, PhiN[1]=654.265063, PhiN[2]=67.775393, State=28;
			//PhiC[1]=416.583354, PhiC[2]=92.872928, PhiN[1]=756.654069, PhiN[2]=530.170519, State=27;
			//PhiC[1]=14.129171, PhiC[2]=936.990153, PhiN[1]=763.860329, PhiN[2]=309.207361, State=26;
			//PhiC[1]=189.087523, PhiC[2]=151.161218, PhiN[1]=496.195400, PhiN[2]=940.222972, State=25;
			//PhiC[1]=930.5060554743, PhiC[2]=384.1185640693, PhiN[1]=654.2650625706, PhiN[2]=67.7753933892, State=24;
			//PhiC[1]=266.8783695698, PhiC[2]=986.6554801464, PhiN[1]=277.8047057986, PhiN[2]=629.9138802290, State=23;
			//PhiC[1]=930.5060554743, PhiC[2]=384.1185640693, PhiN[1]=654.2650625706, PhiN[2]=67.7753933892, State=22;
			//PhiC[1]=416.5833535790, PhiC[2]=92.8729281500, PhiN[1]=756.6540690660, PhiN[2]=530.1705194116, State=21;
			//PhiC[1]=722.9377619028, PhiC[2]=671.4782236814, PhiN[1]=384.0322237611, PhiN[2]=632.0030775070, State=20;
			//PhiC[1]=16.8518330995, PhiC[2]=240.6709003747, PhiN[1]=591.5224537849, PhiN[2]=625.2244701385, State=19;
			//PhiC[1]=73.6131965667, PhiC[2]=273.4372557402, PhiN[1]=897.7586056590, PhiN[2]=275.6319370568, State=18;
			//PhiC[1]=1.5944474395, PhiC[2]=389.4344177544, PhiN[1]=696.5465664864, PhiN[2]=12.3049105685, State=17;
			//PhiC[1]=1.5944474395, PhiC[2]=389.4344177544, PhiN[1]=696.5465664864, PhiN[2]=12.3049105685, State = 16;
			//PhiC[1]=266.8783695698, PhiC[2]=986.6554801464, PhiN[1]=277.8047057986, PhiN[2]=629.9138802290, State = 15;
			//PhiC[1]=416.5833535790, PhiC[2]=92.8729281500, PhiN[1]=756.6540690660, PhiN[2]=530.1705194116, State=14;
			//PhiC[1]=930.5060554743, PhiC[2]=384.1185640693, PhiN[1]=654.2650625706, PhiN[2]=67.7753933892, State=13;
			//PhiC[1]=941.0389776230, PhiC[2]=761.7527322173, PhiN[1]=154.5662565678, PhiN[2]=142.6779553294, State = 12;
			//PhiC[1]=930.5060554743, PhiC[2]=384.1185640693, PhiN[1]=654.2650625706, PhiN[2]=67.7753933892, State=11;
			//PhiC[1]=16.8518330995, PhiC[2]=240.6709003747, PhiN[1]=591.5224537849, PhiN[2]=625.2244701385, State=10;
			//PhiC[1]=29.3468071185, PhiC[2]=587.5994479060, PhiN[1]=565.3337710500, PhiN[2]=471.7910734415, State=9;
			//PhiC[1]=576.5967749357, PhiC[2]=217.1697681397, PhiN[1]=836.7612619996, PhiN[2]=935.3063384295, State=8;
			//PhiC[1]=416.5833535790, PhiC[2]=92.8729281500, PhiN[1]=756.6540690660, PhiN[2]=530.1705194116, State=7;
			//PhiC[1]=416.5833535790, PhiC[2]=92.8729281500, PhiN[1]=756.6540690660, PhiN[2]=530.1705194116, State=6;
			//PhiC[1]=16.8518330995, PhiC[2]=240.6709003747, PhiN[1]=591.5224537849, PhiN[2]=625.2244701385, State=5;
			//PhiC[1]=930.5060554743, PhiC[2]=384.1185640693, PhiN[1]=654.2650625706, PhiN[2]=67.7753933892, State=4;
			//PhiC[1]=930.5060554743, PhiC[2]=384.1185640693, PhiN[1]=654.2650625706, PhiN[2]=67.7753933892, State=3;
			//PhiC[1]=416.5833535790, PhiC[2]=92.8729281500, PhiN[1]=756.6540690660, PhiN[2]=530.1705194116, State=2;
			//PhiC[1]=930.5060554743, PhiC[2]=384.1185640693, PhiN[1]=654.2650625706, PhiN[2]=67.7753933892, State=1;


			double Noise1=10.0;//Noise in the nutrint concentration
			double Noise2=10.0;//Noise in the biomass of bacteria which are PRESENT in the state

			
			for (c=1; c <= n_C; c++)
				C[c] = 0.0;
			for (n=1; n <= n_N; n++)
				N[n] = 0.0;


			
			if(State == 33)		
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass=51.660424496045053 + Noise2,	stSpecies[2].dBiomass=79.083355072994451 + Noise2,
				stSpecies[3].dBiomass=123.698855795127031 + Noise2,	stSpecies[4].dBiomass=69.894543806010333 + Noise2;
			}							


		
			if(State == 32)		
			{	C[1] = (double)(dDiffConst)/stSpecies[1].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;				
				N[1] = 630.440074775209496 + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass=15.680236476604424 + Noise2,	stSpecies[2].dBiomass= 0.0 + Noise2;
				stSpecies[3].dBiomass=55.279539841297606 + Noise2,	stSpecies[4].dBiomass=90.955456117767653 + Noise2;
			}		
			
		
			if(State == 31)		
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[4].LambdaC + Noise1;				
				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = 142.736301530885612 + Noise1;

				stSpecies[1].dBiomass=240.528601860943098 + Noise2,	stSpecies[2].dBiomass= 48.389926660549861 + Noise2;
				stSpecies[3].dBiomass = 0.0 + Noise2,	stSpecies[4].dBiomass=71.926462222194942 + Noise2;
			}	


			if(State == 30)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = 54.539455112677231 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[3].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass=0.0 + Noise2,	stSpecies[2].dBiomass=1.020446338494576 + Noise2,
				stSpecies[3].dBiomass=153.240254176003361 + Noise2,	stSpecies[4].dBiomass=1.239263344534280 + Noise2;
			}			
			
			
			if(State == 29)
			{	C[1] = 152.325656747747900 + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[2].LambdaN + Noise1;

				stSpecies[1].dBiomass=96.107698007910784 + Noise2,	stSpecies[2].dBiomass=52.740182224501055 + Noise2,
				stSpecies[3].dBiomass=23.985932800915659 + Noise2,	stSpecies[4].dBiomass=0.0 + Noise2;
			}			

	
			
			if(State == 28)
			{	C[1] = 919.916150257224217 + Noise1;
				C[2] = 77.866810665455091 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[3].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[2].LambdaN + Noise1;

				stSpecies[1].dBiomass=0.0 + Noise2,	stSpecies[2].dBiomass=6.777539439909746 + Noise2,
				stSpecies[3].dBiomass=143.938322734890249 + Noise2,	stSpecies[4].dBiomass=0.0 + Noise2;
			}		

		
			if(State == 27)
			{	C[1] = (double)(dDiffConst)/stSpecies[1].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[4].LambdaC + Noise1;	

				N[1] = 185.780577174674931 + Noise1;
				N[2] = 508.132875748697984 + Noise1;

				stSpecies[1].dBiomass=154.135842810670624 + Noise2,	stSpecies[2].dBiomass=0.0 + Noise2,
				stSpecies[3].dBiomass=0.0 + Noise2,	stSpecies[4].dBiomass=13.002209996356326 + Noise2;
			}		

			//From here the IS starts
			if(State == 26)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;	

				N[1] = 731.845241307953529 + Noise1;	
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;	

				stSpecies[1].dBiomass=0.0,	stSpecies[2].dBiomass=9.042670396457334 + Noise2,
				stSpecies[3].dBiomass=7.043417811273002 + Noise2,	stSpecies[4].dBiomass=129.080601221175186 + Noise2;
			}
			
			if(State == 25)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;	

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;	
				N[2] = 539.234089407511419 + Noise1;	

				stSpecies[1].dBiomass=46.780226325409856 + Noise2,	stSpecies[2].dBiomass=40.098864745720761 + Noise2,
				stSpecies[3].dBiomass=71.045772329972721 + Noise2,	stSpecies[4].dBiomass=0.0;
			}

		
			if(State == 24)
			{	C[1] = 606.617901722761871 + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass=119.838634212834421 + Noise2,	stSpecies[2].dBiomass=0.0,
				stSpecies[3].dBiomass=46.292032348074144 + Noise2,	stSpecies[4].dBiomass=39.987480322127155 + Noise2;
			}						
			
			if(State == 23)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = 62.397306284005026 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass=75.007273546150145 + Noise2,	stSpecies[2].dBiomass=41.059852001277335 + Noise2,
				stSpecies[3].dBiomass=0.0,	stSpecies[4].dBiomass=129.396072017367430 + Noise2;
			}								

			if(State == 22)
			{	C[1] = 930.506055474281311 + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = 443.846740455798454 + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass=0.0,	stSpecies[2].dBiomass=0.0,
				stSpecies[3].dBiomass=46.292040427530623 + Noise2,	stSpecies[4].dBiomass=39.987480322127155 + Noise2;
			}												

			if(State == 21)
			{	C[1] = 9.216281858962702+ Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] =  530.170519411563873+ Noise1;

				stSpecies[1].dBiomass=150.725814478406619 + Noise2,	stSpecies[2].dBiomass=0.0,
				stSpecies[3].dBiomass=43.650276119785993 + Noise2,	stSpecies[4].dBiomass=0.0;
			}

			if(State == 20)
			{	C[1] =442.698037855294046 + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[4].LambdaC + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = 472.669263922869391 + Noise1;

				stSpecies[1].dBiomass=103.688693090658688 + Noise2,	stSpecies[2].dBiomass=0.0,
				stSpecies[3].dBiomass=0.0,	stSpecies[4].dBiomass=94.006951715635196 + Noise2;
			}			
			
			if(State == 19)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = 77.361925192050535 + Noise1;
				N[2] = 517.372736694637979 + Noise1;

				stSpecies[1].dBiomass=0.0,	stSpecies[2].dBiomass= 10.785172942612174 + Noise2,
				stSpecies[3].dBiomass= 113.115322889183872 + Noise2,	stSpecies[4].dBiomass=0.0;
			}						
			
			if(State == 18)
			{	C[1] = (double)(dDiffConst)/stSpecies[1].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = 212.719895260275734 + Noise1;
				N[2] = 275.631937056779861 + Noise1;

				stSpecies[1].dBiomass= 27.236880886850205 + Noise2,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 128.515509871915270 + Noise2,	stSpecies[4].dBiomass=0.0;
			}									


			if(State == 17)
			{	C[1] = (double)(dDiffConst)/stSpecies[1].LambdaC + Noise1;
				C[2] = 337.578007253521719 + Noise1;

				N[1] = 694.361582929946053 + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass= 0.589945560231405 + Noise2,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 7.259897646147753 + Noise2;
			}
			
			if(State == 16)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = 380.582537304629000 + Noise1;

				N[1] = 696.546566486358643 + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 1.020446338494576 + Noise2,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 1.239263724680173 + Noise2;
			}

			if(State == 15)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = 986.655480146408081 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = 219.315395938166148 + Noise1;

				stSpecies[1].dBiomass= 75.007273546150145 + Noise2,	stSpecies[2].dBiomass= 41.059843722051937 + Noise2,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 0.0;
			}

			if(State == 14)
			{	C[1] = 333.744209920987487 + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = 558.243734641577817 + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[2].LambdaN + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 53.017048781093841 + Noise2,
				stSpecies[3].dBiomass= 43.650276119785993 + Noise2,	stSpecies[4].dBiomass= 0.0;
			}

			if(State == 13)
			{	C[1] = 442.479494569769145 + Noise1;
				C[2] = 384.118564069271088 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[2].LambdaN + Noise1;

				stSpecies[1].dBiomass= 176.651573913557002 + Noise2,	stSpecies[2].dBiomass= 6.777539439909746 + Noise2,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 0.0;
			}

			if(State == 12)
			{	C[1] = 941.038977622985840 + Noise1;
				C[2] = 88.116876047457822 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[3].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 34.004576260666539 + Noise2,	stSpecies[4].dBiomass= 84.179989902478439 + Noise2;
			}

			if(State == 11)
			{	C[1] = 453.069399786826239 + Noise1;
				C[2] = 98.493682695552309 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass= 176.651573913557002 + Noise2,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 39.987480322127155 + Noise2;
			}

			if(State == 10)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[4].LambdaC + Noise1;

				N[1] = 591.522453784942627 + Noise1;
				N[2] = 460.264384657115841 + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 10.785172942612174 + Noise2,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 33.693926195902172 + Noise2;
			}


			if(State == 9)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = 322.975110543097685 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[3].LambdaN + Noise1;
				N[2] = 283.971505084342311 + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 18.781956136033909 + Noise2,
				stSpecies[3].dBiomass= 124.373437381194080 + Noise2,	stSpecies[4].dBiomass= 0.0;
			}
			
			
			if(State == 8)
			{	C[1] = (double)(dDiffConst)/stSpecies[1].LambdaC + Noise1;
				C[2] = 217.169768139719963 + Noise1;

				N[1] = 46.610184529206549 + Noise1;
				N[2] = 935.306338429450989 + Noise1;

				stSpecies[1].dBiomass= 213.340809475644960 + Noise2,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 0.0;
			}
			
			if(State == 7)
			{	C[1] = 333.744209920987487 + Noise1;
				C[2] = 92.872928149998188 + Noise1;

				N[1] = 756.654069066047668 + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[2].LambdaN + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 53.017052731172022 + Noise2,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 0.0;
			}


			if(State == 6)
			{	C[1] = 416.583353579044342 + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[3].LambdaC + Noise1;

				N[1] = 558.243712498946252 + Noise1;
				N[2] = 530.170519411563873 + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 43.650276119785993 + Noise2,	stSpecies[4].dBiomass= 0.0;
			}

			if(State == 5)
			{	C[1] = (double)(dDiffConst)/stSpecies[2].LambdaC + Noise1;
				C[2] = 240.670900374650955 + Noise1;

				N[1] = 591.522453784942627 + Noise1;
				N[2] = 517.372736694637979 + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 10.785172942612174 + Noise2,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 0.0;
			}

			if(State == 4)
			{	C[1] = 930.506055474281311 + Noise1;
				C[2] = 98.493715013378335 + Noise1;

				N[1] = 654.265062570571899 + Noise1;
				N[2] = (double)(dDiffConst)/stSpecies[4].LambdaN + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 39.987480322127155 + Noise2;
			}

			if(State == 3)
			{	C[1] = 930.506055474281311 + Noise1;
				C[2] = 77.866830164073406 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[3].LambdaN + Noise1;
				N[2] = 67.775393389165401 + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 143.938312985581092 + Noise2,	stSpecies[4].dBiomass= 0.0;
			}

			if(State == 2)
			{	C[1] = 416.583353579044342 + Noise1;
				C[2] = (double)(dDiffConst)/stSpecies[4].LambdaC + Noise1;

				N[1] = 756.654069066047668 + Noise1;
				N[2] = 508.132874364783504 + Noise1;

				stSpecies[1].dBiomass= 0.0,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 13.002209996356326 + Noise2;
			}

			if(State == 1)
			{	C[1] = 453.069360789589609 + Noise1;
				C[2] = 384.118564069271088 + Noise1;

				N[1] = (double)(dDiffConst)/stSpecies[1].LambdaN + Noise1;
				N[2] = 67.775393389165401 + Noise1;

				stSpecies[1].dBiomass= 176.651573913557002 + Noise2,	stSpecies[2].dBiomass= 0.0,
				stSpecies[3].dBiomass= 0.0,	stSpecies[4].dBiomass= 0.0;
			}



					printf("\n Initially:\n");
					for( b=1; b<= n_C*n_N; b++)
						printf("B[%ld]=%f\n", b, stSpecies[b].dBiomass);
					printf("C1=%f, C2=%f, N1=%f, N2=%f\n", C[1], C[2], N[1], N[2]);
					
						

					vDyneig();	//integrate the system

					printf("\n Finally:\n");
					for( b=1; b<= n_C*n_N; b++)
						printf("B[%ld]=%f\n", b, stSpecies[b].dBiomass);
					printf("C1=%f, C2=%f, N1=%f, N2=%f\n", C[1], C[2], N[1], N[2]);
					
					printf("\nState %ld now is:	", State);

					//Checking the sign of the absent species
					for(b=1; b <= n_C*n_N; b++)
					{	
						if(stSpecies[b].dBiomass < 1e-7)
		 				{	stSpecies[b].dBiomass = 0.0;
		
							if ( (stSpecies[b].LambdaC * C[G[b][1]]) < (stSpecies[b].LambdaN * N[G[b][2]]) )
								Min = (stSpecies[b].LambdaC * C[G[b][1]]);
							else
								Min = (stSpecies[b].LambdaN * N[G[b][2]]);
							
							if ( (Min - dDiffConst) > 0.0)
				  				NoOfPlus += 1;						
							printf("0  ");				  				
		  				}
		  				
						else
		 				{		
							if ( (stSpecies[b].LambdaC * C[G[b][1]]) < (stSpecies[b].LambdaN * N[G[b][2]]) )
								printf("-1 ");
							else
								printf("1  ");
							
							if ( (Min - dDiffConst) > 0.0)
				  				NoOfPlus += 1;						
		  				}		  				
		  				
					}
					
					printf("\nNoOfPlus = %ld\n", NoOfPlus);

		fclose(fFinalStats);

	//printf("\nMinimum Resource Mass (throughout the experiment) = %.5f,	", dMinResourceMass);
	return 0;
}
