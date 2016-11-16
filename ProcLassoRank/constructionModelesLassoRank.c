#include "EMGrank.h"
#include "constructionModelesLassoRank.h"
#include <gsl/gsl_linalg.h>
#include <omp.h>
#include "omp_num_threads.h"

// TODO: comment on constructionModelesLassoRank purpose
void constructionModelesLassoRank(
	// IN parameters 
	const Real* Pi,    // parametre initial des proportions
	const Real* Rho,   // parametre initial de variance renormalisé
	Int mini,       // nombre minimal d'itérations dans l'algorithme EM	
	Int maxi,       // nombre maximal d'itérations dans l'algorithme EM
	const Real* X,      // régresseurs
	const Real* Y,      // réponse
	Real tau,     // seuil pour accepter la convergence
	const Int* A1,     // matrice des coefficients des parametres selectionnes
	Int rangmin, 	//rang minimum autorisé
	Int rangmax,	//rang maximum autorisé
	// OUT parameters (all pointers, to be modified)
	Real* phi,        // estimateur ainsi calculé par le Lasso
    Real* lvraisemblance,      // estimateur ainsi calculé par le Lasso
    // additional size parameters
	mwSize n,              // taille de l'echantillon                
	mwSize p,              // nombre de covariables
	mwSize m,              // taille de Y (multivarié)
	mwSize k,              // nombre de composantes
	mwSize L)              // taille de glambda
{
	//On cherche les rangs possiblement intéressants
	Int deltaRank = rangmax-rangmin+1;
	mwSize Size = (mwSize)pow(deltaRank,k);
	Int* Rank = (Int*)malloc(Size*k*sizeof(Int));
    for (mwSize r=0; r<k; r++)
    {
		//On veut le tableau de toutes les combinaisons de rangs possibles
		//Dans la première colonne : on répète (rangmax-rangmin)^(k-1) chaque chiffre : ca remplit la colonne
		//Dans la deuxieme : on répète (rangmax-rangmin)^(k-2) chaque chiffre, et on fait ca (rangmax-rangmin)^2 fois 
		//...
		//Dans la dernière, on répète chaque chiffre une fois, et on fait ca (rangmin-rangmax)^(k-1) fois.
		Int indexInRank = 0;
		Int value = 0;
		while (indexInRank < Size)
		{
			for (Int u=0; u<pow(deltaRank,k-r-1); u++)
				Rank[(indexInRank++)*k+r] = rangmin + value;
			value = (value+1) % deltaRank;
		}
	}
	
	//Initialize phi to zero, because unactive variables won't be assigned
	for (mwSize i=0; i<p*m*k*L*Size; i++)
		phi[i] = 0.0;
	
	//initiate parallel section
	mwSize lambdaIndex;
	omp_set_num_threads(OMP_NUM_THREADS);
	#pragma omp parallel default(shared) private(lambdaIndex)
	{
	#pragma omp for schedule(dynamic,CHUNK_SIZE) nowait
	for (lambdaIndex=0; lambdaIndex<L; lambdaIndex++)
	{
		//On ne garde que les colonnes actives : active sera l'ensemble des variables informatives
		Int* active = (Int*)malloc(p*sizeof(Int));
		mwSize longueurActive = 0;
		for (Int j=0; j<p; j++)
		{
			if (A1[j*L+lambdaIndex] != 0)
				active[longueurActive++] = A1[j*L+lambdaIndex] - 1;
		}
		
		if (longueurActive == 0)
			continue;
		
		//from now on, longueurActive > 0
		Real* phiLambda = (Real*)malloc(longueurActive*m*k*sizeof(Real));
		Real LLF;
		for (Int j=0; j<Size; j++)
		{
			//[phiLambda,LLF] = EMGrank(Pi(:,lambdaIndex),Rho(:,:,:,lambdaIndex),mini,maxi,X(:,active),Y,tau,Rank(j,:));
			Int* rank = (Int*)malloc(k*sizeof(Int));
			for (mwSize r=0; r<k; r++)
				rank[r] = Rank[j*k+r];
			Real* Xactive = (Real*)malloc(n*longueurActive*sizeof(Real));
			for (mwSize i=0; i<n; i++)
			{
				for (mwSize jj=0; jj<longueurActive; jj++)
					Xactive[i*longueurActive+jj] = X[i*p+active[jj]];
			}
			Real* PiLambda = (Real*)malloc(k*sizeof(Real));
			for (mwSize r=0; r<k; r++)
				PiLambda[r] = Pi[r*L+lambdaIndex];
			Real* RhoLambda = (Real*)malloc(m*m*k*sizeof(Real));
			for (mwSize u=0; u<m; u++)
			{
				for (mwSize v=0; v<m; v++)
				{
					for (mwSize r=0; r<k; r++)
						RhoLambda[u*m*k+v*k+r] = Rho[u*m*k*L+v*k*L+r*L+lambdaIndex];
				}
			}
			EMGrank(PiLambda,RhoLambda,mini,maxi,Xactive,Y,tau,rank,
				phiLambda,&LLF,
				n,longueurActive,m,k);
			free(rank);
			free(Xactive);
			free(PiLambda);
			free(RhoLambda);
			//lvraisemblance((lambdaIndex-1)*Size+j,:) = [LLF, dot(Rank(j,:), length(active)-Rank(j,:)+m)];
			lvraisemblance[(lambdaIndex*Size+j)*2] = LLF;
			//dot(Rank(j,:), length(active)-Rank(j,:)+m)
			Real dotProduct = 0.0;
			for (mwSize r=0; r<k; r++)
				dotProduct += Rank[j*k+r] * (longueurActive-Rank[j*k+r]+m);
			lvraisemblance[(lambdaIndex*Size+j)*2+1] = dotProduct;
			//phi(active,:,:,(lambdaIndex-1)*Size+j) = phiLambda;
			for (mwSize jj=0; jj<longueurActive; jj++)
			{
				for (mwSize mm=0; mm<m; mm++)
				{
					for (mwSize r=0; r<k; r++)
						phi[active[jj]*m*k*L*Size+mm*k*L*Size+r*L*Size+(lambdaIndex*Size+j)] = phiLambda[jj*m*k+mm*k+r];
				}
			}
		}
		free(active);
		free(phiLambda);
	}
	}
	free(Rank);
}
