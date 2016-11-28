#include "selectiontotale.h"
#include "EMGLLF.h"
#include <omp.h>
#include "omp_num_threads.h"

// Main job on raw inputs (after transformation from mxArray)
void selectiontotale(
	// IN parameters 
	const Real* phiInit,   // parametre initial de moyenne renormalisé
	const Real* rhoInit,   // parametre initial de variance renormalisé
	const Real* piInit,    // parametre initial des proportions
	const Real* gamInit,   // paramètre initial des probabilités a posteriori de chaque échantillon
	Int mini,       // nombre minimal d'itérations dans lambdaIndex'algorithme EM	
	Int maxi,       // nombre maximal d'itérations dans lambdaIndex'algorithme EM
	Real gamma,   // valeur de gamma : puissance des proportions dans la pénalisation pour un Lasso adaptatif
	const Real* glambda, // valeur des paramètres de régularisation du Lasso
	const Real* X,      // régresseurs
	const Real* Y,      // réponse
	Real seuil,   // seuil pour prendre en compte une variable
	Real tau,     // seuil pour accepter la convergence
	// OUT parameters (all pointers, to be modified)
	Int* A1,         // matrice des coefficients des parametres selectionnes
	Int* A2,         // matrice des coefficients des parametres non selectionnes
	Real* Rho,        // estimateur ainsi calculé par le Lasso
	Real* Pi,          // estimateur ainsi calculé par le Lasso
	// additional size parameters
	mwSize n,              // taille de lambdaIndex'echantillon                
	mwSize p,              // nombre de covariables
	mwSize m,              // taille de Y (multivarié)
	mwSize k,              // nombre de composantes
	mwSize L)             // taille de glambda
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
	mwSize lambdaIndex;
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
		EMGLLF(phiInit,rhoInit,piInit,gamInit,mini,maxi,gamma,glambda[lambdaIndex],X,Y,tau,
			phi,rho,pi,LLF,S,
			n,p,m,k);
		free(LLF);
		free(S);
		
		//Si un des coefficients est supérieur au seuil, on garde cette variable
		mwSize* selectedVariables = (mwSize*)calloc(p*m,sizeof(mwSize));
		mwSize* discardedVariables = (mwSize*)calloc(p*m,sizeof(mwSize));
		int atLeastOneSelectedVariable = 0;
		for (mwSize j=0; j<p; j++)
		{
			mwSize cpt = 0;
			mwSize cpt2 = 0;
			for (mwSize jj=0; jj<m; jj++)
			{
				Real maxPhi = 0.0;
				for (mwSize r=0; r<k; r++)
				{
					if (fabs(phi[j*m*k+jj*k+r]) > maxPhi)
						maxPhi = fabs(phi[j*m*k+jj*k+r]);
				}
				if (maxPhi > seuil)
				{
					selectedVariables[j*m+cpt] = jj+1;
					atLeastOneSelectedVariable = 1;
					cpt++;
				}
				else
				{
					discardedVariables[j*m+cpt2] = jj+1;
					cpt2++;
				}
			}
		}
		free(phi);
		
		if (atLeastOneSelectedVariable)
		{
			Int* vec = (Int*)malloc(p*sizeof(Int));
			mwSize vecSize = 0;
			for (mwSize j=0; j<p; j++)
			{
				if (selectedVariables[j*m+0] != 0)
					vec[vecSize++] = j;
			}
			
			//A1
			for (mwSize j=0; j<p; j++)
				A1[j*(m+1)*L+0*L+lambdaIndex] = (j < vecSize ? vec[j]+1 : 0);
			for (mwSize j=0; j<vecSize; j++)
			{
				for (mwSize jj=1; jj<=m; jj++)
					A1[j*(m+1)*L+jj*L+lambdaIndex] = selectedVariables[vec[j]*m+jj-1];
			}
			//A2
			for (mwSize j=0; j<p; j++)
				A2[j*(m+1)*L+0*L+lambdaIndex] = j+1;
			for (mwSize j=0; j<p; j++)
			{
				for (mwSize jj=1; jj<=m; jj++)
					A2[j*(m+1)*L+jj*L+lambdaIndex] = discardedVariables[j*m+jj-1];
			}
			//Rho
			for (mwSize j=0; j<m; j++)
			{
				for (mwSize jj=0; jj<m; jj++)
				{
					for (mwSize r=0; r<k; r++)
						Rho[j*m*k*L+jj*k*L+r*L+lambdaIndex] = rho[j*m*k+jj*k+r];
				}
			}
			//Pi
			for (mwSize r=0; r<k; r++)
				Pi[r*L+lambdaIndex] = pi[r];
			free(vec);
		}
		
		free(selectedVariables);
		free(discardedVariables);
		free(rho);
		free(pi);
	}
	}
}
