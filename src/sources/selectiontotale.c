#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "EMGLLF.h"
#include "utils.h"

// Main job on raw inputs (after transformation from mxArray)
void selectiontotale_core(
	// IN parameters
	const Real* phiInit, // parametre initial de moyenne renormalisé
	const Real* rhoInit, // parametre initial de variance renormalisé
	const Real* piInit,// parametre initial des proportions
	const Real* gamInit, // paramètre initial des probabilités a posteriori de chaque échantillon
	int mini, // nombre minimal d'itérations dans lambdaIndex'algorithme EM
	int maxi, // nombre maximal d'itérations dans lambdaIndex'algorithme EM
	Real gamma, // valeur de gamma : puissance des proportions dans la pénalisation pour un Lasso adaptatif
	const Real* glambda, // valeur des paramètres de régularisation du Lasso
	const Real* X,// régresseurs
	const Real* Y,// réponse
	Real seuil, // seuil pour prendre en compte une variable
	Real tau, // seuil pour accepter la convergence
	// OUT parameters (all pointers, to be modified)
	int* A1, // matrice des coefficients des parametres selectionnes
	int* A2, // matrice des coefficients des parametres non selectionnes
	Real* Rho,// estimateur ainsi calculé par le Lasso
	Real* Pi,// estimateur ainsi calculé par le Lasso
	// additional size parameters
	int n,// taille de lambdaIndex'echantillon
	int p,// nombre de covariables
	int m,// taille de Y (multivarié)
	int k,// nombre de composantes
	int L) // taille de glambda
{
	// Fill outputs with zeros: they might not be assigned
	for (int u=0; u<p*(m+1)*L; u++)
	{
		A1[u] = 0;
		A2[u] = 0;
	}
	for (int u=0; u<m*m*k*L; u++)
		Rho[u] = 0.0;
	for (int u=0; u<k*L; u++)
		Pi[u] = 0.0;

	//initiate parallel section
	int lambdaIndex;
	omp_set_num_threads(OMP_NUM_THREADS);
	#pragma omp parallel default(shared) private(lambdaIndex)
	{
	#pragma omp for schedule(dynamic,CHUNK_SIZE) nowait
	for (lambdaIndex=0; lambdaIndex<L; lambdaIndex++)
	{
		//allocate output variables
		Real* phi = (Real*)malloc(p*m*k*sizeof(Real));
		Real* rho = (Real*)malloc(m*m*k*sizeof(Real));
		Real* pi = (Real*)malloc(k*sizeof(Real));
		Real* LLF = (Real*)malloc(maxi*sizeof(Real));
		Real* S = (Real*)malloc(p*m*k*sizeof(Real));
		EMGLLF_core(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda[lambdaIndex],X,Y,tau,
			phi,rho,pi,LLF,S,
			n,p,m,k);
		free(LLF);
		free(S);

		//Si un des coefficients est supérieur au seuil, on garde cette variable
		int* selectedVariables = (int*)calloc(p*m,sizeof(int));
		int* discardedVariables = (int*)calloc(p*m,sizeof(int));
		int atLeastOneSelectedVariable = 0;
		for (int j=0; j<p; j++)
		{
			int cpt = 0;
			int cpt2 = 0;
			for (int jj=0; jj<m; jj++)
			{
				Real maxPhi = 0.0;
				for (int r=0; r<k; r++)
				{
					if (fabs(phi[ai(j,jj,r,p,m,k)]) > maxPhi)
						maxPhi = fabs(phi[ai(j,jj,r,p,m,k)]);
				}
				if (maxPhi > seuil)
				{
					selectedVariables[mi(j,cpt,p,m)] = jj+1;
					atLeastOneSelectedVariable = 1;
					cpt++;
				}
				else
				{
					discardedVariables[mi(j,cpt2,p,m)] = jj+1;
					cpt2++;
				}
			}
		}
		free(phi);

		if (atLeastOneSelectedVariable)
		{
			int* vec = (int*)malloc(p*sizeof(int));
			int vecSize = 0;
			for (int j=0; j<p; j++)
			{
				if (selectedVariables[mi(j,0,p,m)] != 0)
					vec[vecSize++] = j;
			}
			
			//A1
			for (int j=0; j<p; j++)
				A1[ai(j,0,lambdaIndex,p,m+1,L)] = (j < vecSize ? vec[j]+1 : 0);
			for (int j=0; j<vecSize; j++)
			{
				for (int jj=1; jj<=m; jj++)
					A1[ai(j,jj,lambdaIndex,p,m+1,L)] = selectedVariables[mi(vec[j],jj-1,p,m)];
			}
			//A2
			for (int j=0; j<p; j++)
				A2[ai(j,0,lambdaIndex,p,m+1,L)] = j+1;
			for (int j=0; j<p; j++)
			{
				for (int jj=1; jj<=m; jj++)
					A2[ai(j,jj,lambdaIndex,p,m+1,L)] = discardedVariables[mi(j,jj-1,p,m)];
			}
			//Rho
			for (int j=0; j<m; j++)
			{
				for (int jj=0; jj<m; jj++)
				{
					for (int r=0; r<k; r++)
						Rho[ai4(j,jj,r,lambdaIndex,m,m,k,L)] = rho[ai(j,jj,r,m,m,k)];
				}
			}
			//Pi
			for (int r=0; r<k; r++)
				Pi[mi(r,lambdaIndex,k,L)] = pi[r];
			free(vec);
		}

		free(selectedVariables);
		free(discardedVariables);
		free(rho);
		free(pi);
	}
	}
}
