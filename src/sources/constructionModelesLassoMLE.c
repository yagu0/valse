#include "EMGLLF.h"
#include "utils.h"
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <omp.h>

// TODO: comment on constructionModelesLassoMLE purpose
void constructionModelesLassoMLE(
	// IN parameters
	const double* phiInit, // parametre initial de moyenne renormalisé
	const double* rhoInit, // parametre initial de variance renormalisé
	const double* piInit,// parametre initial des proportions
	const double* gamInit, // paramètre initial des probabilités a posteriori de chaque échantillon
	int mini,// nombre minimal d'itérations dans l'algorithme EM
	int maxi,// nombre maximal d'itérations dans l'algorithme EM
	double gamma,// valeur de gamma : puissance des proportions dans la pénalisation pour un Lasso adaptatif
	const double* glambda, // valeur des paramètres de régularisation du Lasso
	const double* X, // régresseurs
	const double* Y, // réponse
	double seuil,// seuil pour prendre en compte une variable
	double tau,// seuil pour accepter la convergence
	const int* A1, // matrice des coefficients des parametres selectionnes
	const int* A2, // matrice des coefficients des parametres non selectionnes
	// OUT parameters
	double* phi,// estimateur ainsi calculé par le Lasso
	double* rho,// estimateur ainsi calculé par le Lasso
	double* pi, // estimateur ainsi calculé par le Lasso
	double* lvraisemblance, // estimateur ainsi calculé par le Lasso
	// additional size parameters
	int n, // taille de l'echantillon
	int p, // nombre de covariables
	int m, // taille de Y (multivarié)
	int k, // nombre de composantes
	int L) // taille de glambda
{
	//preparation: phi = 0
	for (int u=0; u<p*m*k*L; u++)
		phi[u] = 0.0;

	//initiate parallel section
	int lambdaIndex;
	omp_set_num_threads(OMP_NUM_THREADS);
	#pragma omp parallel default(shared) private(lambdaIndex)
	{
	#pragma omp for schedule(dynamic,CHUNK_SIZE) nowait
	for (lambdaIndex=0; lambdaIndex<L; lambdaIndex++)
	{
		//~ a = A1(:,1,lambdaIndex);
		//~ a(a==0) = [];
		int* a = (int*)malloc(p*sizeof(int));
		int lengthA = 0;
		for (int j=0; j<p; j++)
		{
			if (A1[ai(j,0,lambdaIndex,p,m+1,L)] != 0)
				a[lengthA++] = A1[ai(j,0,lambdaIndex,p,m+1,L)] - 1;
		}
		if (lengthA == 0)
			continue;

		//Xa = X(:,a)
		double* Xa = (double*)malloc(n*lengthA*sizeof(double));
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<lengthA; j++)
				Xa[mi(i,j,n,lengthA)] = X[mi(i,a[j],n,p)];
		}

		//phia = phiInit(a,:,:)
		double* phia = (double*)malloc(lengthA*m*k*sizeof(double));
		for (int j=0; j<lengthA; j++)
		{
			for (int mm=0; mm<m; mm++)
			{
				for (int r=0; r<k; r++)
					phia[ai(j,mm,r,lengthA,m,k)] = phiInit[ai(a[j],mm,r,p,m,k)];
			}
		}

		//[phiLambda,rhoLambda,piLambda,~,~] = EMGLLF(...
		//	phiInit(a,:,:),rhoInit,piInit,gamInit,mini,maxi,gamma,0,X(:,a),Y,tau);
		double* phiLambda = (double*)malloc(lengthA*m*k*sizeof(double));
		double* rhoLambda = (double*)malloc(m*m*k*sizeof(double));
		double* piLambda = (double*)malloc(k*sizeof(double));
		double* LLF = (double*)malloc((maxi+1)*sizeof(double));
		double* S = (double*)malloc(lengthA*m*k*sizeof(double));
		EMGLLF(phia,rhoInit,piInit,gamInit,mini,maxi,gamma,0.0,Xa,Y,tau,
			phiLambda,rhoLambda,piLambda,LLF,S,
			n,lengthA,m,k);
		free(Xa);
		free(phia);
		free(LLF);
		free(S);

		//~ for j=1:length(a)
			//~ phi(a(j),:,:,lambdaIndex) = phiLambda(j,:,:);
		//~ end
		for (int j=0; j<lengthA; j++)
		{
			for (int mm=0; mm<m; mm++)
			{
				for (int r=0; r<k; r++)
					phi[ai4(a[j],mm,r,lambdaIndex,p,m,k,L)] = phiLambda[ai(j,mm,r,p,m,k)];
			}
		}
		free(phiLambda);
		//~ rho(:,:,:,lambdaIndex) = rhoLambda;
		for (int u=0; u<m; u++)
		{
			for (int v=0; v<m; v++)
			{
				for (int r=0; r<k; r++)
					rho[ai4(u,v,r,lambdaIndex,m,m,k,L)] = rhoLambda[ai(u,v,r,m,m,k)];
			}
		}
		free(rhoLambda);
		//~ pi(:,lambdaIndex) = piLambda;
		for (int r=0; r<k; r++)
			pi[mi(r,lambdaIndex,k,L)] = piLambda[r];
		free(piLambda);

		int dimension = 0;
		int* b = (int*)malloc(m*sizeof(int));
		for (int j=0; j<p; j++)
		{
			//~ b = A2(j,2:end,lambdaIndex);
			//~ b(b==0) = [];
			int lengthB = 0;
			for (int mm=0; mm<m; mm++)
			{
				if (A2[ai(j,mm+1,lambdaIndex,p,m+1,L)] != 0)
					b[lengthB++] = A2[ai(j,mm+1,lambdaIndex,p,m+1,L)] - 1;
			}
			//~ if length(b) > 0
				//~ phi(A2(j,1,lambdaIndex),b,:,lambdaIndex) = 0.0;
			//~ end
			if (lengthB > 0)
			{
				for (int mm=0; mm<lengthB; mm++)
				{
					for (int r=0; r<k; r++)
						phi[ai( A2[ai4(j,0,lambdaIndex,p,m+1,L)]-1, b[mm], r, lambdaIndex, p, m, k, L)] = 0.0;
				}
			}

			//~ c = A1(j,2:end,lambdaIndex);
			//~ c(c==0) = [];
			//~ dimension = dimension + length(c);
			for (int mm=0; mm<m; mm++)
			{
				if (A1[ai(j,mm+1,lambdaIndex,p,m+1,L)] != 0)
					dimension++;
			}
		}
		free(b);

		int signum;
		double* densite = (double*)calloc(L*n,sizeof(double));
		double sumLogDensit = 0.0;
		gsl_matrix* matrix = gsl_matrix_alloc(m, m);
		gsl_permutation* permutation = gsl_permutation_alloc(m);
		double* YiRhoR = (double*)malloc(m*sizeof(double));
		double* XiPhiR = (double*)malloc(m*sizeof(double));
		for (int i=0; i<n; i++)
		{
			//~ for r=1:k
				//~ delta = Y(i,:)*rho(:,:,r,lambdaIndex) - (X(i,a)*(phi(a,:,r,lambdaIndex)));
				//~ densite(i,lambdaIndex) = densite(i,lambdaIndex) +...
					//~ pi(r,lambdaIndex)*det(rho(:,:,r,lambdaIndex))/(sqrt(2*PI))^m*exp(-dot(delta,delta)/2.0);
			//~ end
			for (int r=0; r<k; r++)
			{
				//compute det(rho(:,:,r,lambdaIndex)) [TODO: avoid re-computations]
				for (int u=0; u<m; u++)
				{
					for (int v=0; v<m; v++)
						matrix->data[u*m+v] = rho[ai4(u,v,r,lambdaIndex,m,m,k,L)];
				}
				gsl_linalg_LU_decomp(matrix, permutation, &signum);
				double detRhoR = gsl_linalg_LU_det(matrix, signum);

				//compute Y(i,:)*rho(:,:,r,lambdaIndex)
				for (int u=0; u<m; u++)
				{
					YiRhoR[u] = 0.0;
					for (int v=0; v<m; v++)
						YiRhoR[u] += Y[mi(i,v,n,m)] * rho[ai4(v,u,r,lambdaIndex,m,m,k,L)];
				}

				//compute X(i,a)*phi(a,:,r,lambdaIndex)
				for (int u=0; u<m; u++)
				{
					XiPhiR[u] = 0.0;
					for (int v=0; v<lengthA; v++)
						XiPhiR[u] += X[mi(i,a[v],n,p)] * phi[ai4(a[v],u,r,lambdaIndex,p,m,k,L)];
				}
				// On peut remplacer X par Xa dans ce dernier calcul, mais je ne sais pas si c'est intéressant ...

				// compute dotProduct < delta . delta >
				double dotProduct = 0.0;
				for (int u=0; u<m; u++)
					dotProduct += (YiRhoR[u]-XiPhiR[u]) * (YiRhoR[u]-XiPhiR[u]);

				densite[mi(lambdaIndex,i,L,n)] += (pi[mi(r,lambdaIndex,k,L)]*detRhoR/pow(sqrt(2.0*M_PI),m))*exp(-dotProduct/2.0);
			}
			sumLogDensit += log(densite[lambdaIndex*n+i]);
		}
		lvraisemblance[mi(lambdaIndex,0,L,2)] = sumLogDensit;
		lvraisemblance[mi(lambdaIndex,1,L,2)] = (dimension+m+1)*k-1;

		free(a);
		free(YiRhoR);
		free(XiPhiR);
		free(densite);
		gsl_matrix_free(matrix);
		gsl_permutation_free(permutation);
	}
	}
}
