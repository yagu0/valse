#include <stdlib.h>
#include <omp.h>
#include <gsl/gsl_linalg.h>
#include "EMGrank.h"
#include "utils.h"

// TODO: comment on constructionModelesLassoRank purpose
void constructionModelesLassoRank(
	// IN parameters
	const double* Pi,// parametre initial des proportions
	const double* Rho, // parametre initial de variance renormalisé
	int mini, // nombre minimal d'itérations dans l'algorithme EM
	int maxi, // nombre maximal d'itérations dans l'algorithme EM
	const double* X,// régresseurs
	const double* Y,// réponse
	double tau, // seuil pour accepter la convergence
	const int* A1, // matrice des coefficients des parametres selectionnes
	int rangmin,	//rang minimum autorisé
	int rangmax,	//rang maximum autorisé
	// OUT parameters (all pointers, to be modified)
	double* phi,// estimateur ainsi calculé par le Lasso
	double* lvraisemblance,// estimateur ainsi calculé par le Lasso
	// additional size parameters
	int n,// taille de l'echantillon
	int p,// nombre de covariables
	int m,// taille de Y (multivarié)
	int k,// nombre de composantes
	int L)// taille de glambda
{
	//On cherche les rangs possiblement intéressants
	int deltaRank = rangmax-rangmin+1;
	int Size = (int)pow(deltaRank,k);
	int* Rank = (int*)malloc(Size*k*sizeof(int));
for (int r=0; r<k; r++)
{
		//On veut le tableau de toutes les combinaisons de rangs possibles
		//Dans la première colonne : on répète (rangmax-rangmin)^(k-1) chaque chiffre : ca remplit la colonne
		//Dans la deuxieme : on répète (rangmax-rangmin)^(k-2) chaque chiffre, et on fait ca (rangmax-rangmin)^2 fois 
		//...
		//Dans la dernière, on répète chaque chiffre une fois, et on fait ca (rangmin-rangmax)^(k-1) fois.
		int indexInRank = 0;
		int value = 0;
		while (indexInRank < Size)
		{
			for (int u=0; u<pow(deltaRank,k-r-1); u++)
				Rank[mi(indexInRank++,r,Size,k)] = rangmin + value;
			value = (value+1) % deltaRank;
		}
	}

	//Initialize phi to zero, because unactive variables won't be assigned
	for (int i=0; i<p*m*k*L*Size; i++)
		phi[i] = 0.0;

	//initiate parallel section
	int lambdaIndex;
	omp_set_num_threads(OMP_NUM_THREADS);
	#pragma omp parallel default(shared) private(lambdaIndex)
	{
	#pragma omp for schedule(dynamic,CHUNK_SIZE) nowait
	for (lambdaIndex=0; lambdaIndex<L; lambdaIndex++)
	{
		//On ne garde que les colonnes actives : active sera l'ensemble des variables informatives
		int* active = (int*)malloc(p*sizeof(int));
		int longueurActive = 0;
		for (int j=0; j<p; j++)
		{
			if (A1[mi(j,lambdaIndex,p,L)] != 0)
				active[longueurActive++] = A1[mi(j,lambdaIndex,p,L)] - 1;
		}

		if (longueurActive == 0)
			continue;

		//from now on, longueurActive > 0
		double* phiLambda = (double*)malloc(longueurActive*m*k*sizeof(double));
		double LLF;
		for (int j=0; j<Size; j++)
		{
			//[phiLambda,LLF] = EMGrank(Pi(:,lambdaIndex),Rho(:,:,:,lambdaIndex),mini,maxi,X(:,active),Y,tau,Rank(j,:));
			int* rank = (int*)malloc(k*sizeof(int));
			for (int r=0; r<k; r++)
				rank[r] = Rank[mi(j,r,Size,k)];
			double* Xactive = (double*)malloc(n*longueurActive*sizeof(double));
			for (int i=0; i<n; i++)
			{
				for (int jj=0; jj<longueurActive; jj++)
					Xactive[mi(i,jj,n,longueurActive)] = X[mi(i,active[jj],n,p)];
			}
			double* PiLambda = (double*)malloc(k*sizeof(double));
			for (int r=0; r<k; r++)
				PiLambda[r] = Pi[mi(r,lambdaIndex,k,L)];
			double* RhoLambda = (double*)malloc(m*m*k*sizeof(double));
			for (int u=0; u<m; u++)
			{
				for (int v=0; v<m; v++)
				{
					for (int r=0; r<k; r++)
						RhoLambda[ai(uu,v,r,m,m,k)] = Rho[ai4(u,v,r,lambdaIndex,m,m,k,L)];
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
			lvraisemblance[mi(lambdaIndex*Size+j,0,L*Size,2)] = LLF;
			//dot(Rank(j,:), length(active)-Rank(j,:)+m)
			double dotProduct = 0.0;
			for (int r=0; r<k; r++)
				dotProduct += Rank[mi(j,r,Size,k)] * (longueurActive-Rank[mi(j,r,Size,k)]+m);
			lvraisemblance[mi(lambdaIndex*Size+j,1,Size*L,2)] = dotProduct;
			//phi(active,:,:,(lambdaIndex-1)*Size+j) = phiLambda;
			for (int jj=0; jj<longueurActive; jj++)
			{
				for (int mm=0; mm<m; mm++)
				{
					for (int r=0; r<k; r++)
						phi[ai5(active[jj],mm,r,lambdaIndex,j,p,m,k,L,Size)] = phiLambda[jj*m*k+mm*k+r];
				}
			}
		}
		free(active);
		free(phiLambda);
	}
	}
	free(Rank);
}
