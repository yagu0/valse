#include "EMGLLF.h"
#include "utils.h"
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <omp.h>

// TODO: comment on constructionModelesLassoMLE purpose
void constructionModelesLassoMLE_core(
	// IN parameters
	const Real* phiInit, // parametre initial de moyenne renormalisé
	const Real* rhoInit, // parametre initial de variance renormalisé
	const Real* piInit,// parametre initial des proportions
	const Real* gamInit, // paramètre initial des probabilités a posteriori de chaque échantillon
	int mini,// nombre minimal d'itérations dans l'algorithme EM
	int maxi,// nombre maximal d'itérations dans l'algorithme EM
	Real gamma,// valeur de gamma : puissance des proportions dans la pénalisation
	           //pour un Lasso adaptatif
	const Real* glambda, // valeur des paramètres de régularisation du Lasso
	const Real* X, // régresseurs
	const Real* Y, // réponse
	Real seuil,// seuil pour prendre en compte une variable
	Real tau,// seuil pour accepter la convergence
	const int* A1, // matrice des coefficients des parametres selectionnes
	const int* A2, // matrice des coefficients des parametres non selectionnes
	// OUT parameters
	Real* phi,// estimateur ainsi calculé par le Lasso
	Real* rho,// estimateur ainsi calculé par le Lasso
	Real* pi, // estimateur ainsi calculé par le Lasso
	Real* llh, // estimateur ainsi calculé par le Lasso
	// additional size parameters
	int n, // taille de l'echantillon
	int p, // nombre de covariables
	int m, // taille de Y (multivarié)
	int k, // nombre de composantes
	int L) // taille de glambda
{
	//preparation: phi,rho,pi = 0, llh=+Inf
	for (int u=0; u<p*m*k*L; u++)
		phi[u] = 0.;
	for (int u=0; u<m*m*k*L; u++)
		rho[u] = 0.;
	for (int u=0; u<k*L; u++)
		pi[u] = 0.;
	for (int u=0; u<L*2; u++)
		llh[u] = INFINITY;

	//initiate parallel section
	int lambdaIndex;
	omp_set_num_threads(OMP_NUM_THREADS);
	#pragma omp parallel default(shared) private(lambdaIndex)
	{
	#pragma omp for schedule(dynamic,CHUNK_SIZE) nowait
	for (lambdaIndex=0; lambdaIndex<L; lambdaIndex++)
	{
		//a = A1[,1,lambdaIndex] ; a = a[a!=0]
		int* a = (int*)malloc(p*sizeof(int));
		int lengthA = 0;
		for (int j=0; j<p; j++)
		{
			if (A1[ai(j,0,lambdaIndex,p,m+1,L)] != 0)
				a[lengthA++] = A1[ai(j,0,lambdaIndex,p,m+1,L)] - 1;
		}
		if (lengthA == 0)
		{
			free(a);
			continue;
		}

		//Xa = X[,a]
		Real* Xa = (Real*)malloc(n*lengthA*sizeof(Real));
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<lengthA; j++)
				Xa[mi(i,j,n,lengthA)] = X[mi(i,a[j],n,p)];
		}

		//phia = phiInit[a,,]
		Real* phia = (Real*)malloc(lengthA*m*k*sizeof(Real));
		for (int j=0; j<lengthA; j++)
		{
			for (int mm=0; mm<m; mm++)
			{
				for (int r=0; r<k; r++)
					phia[ai(j,mm,r,lengthA,m,k)] = phiInit[ai(a[j],mm,r,p,m,k)];
			}
		}

		//Call to EMGLLF
		Real* phiLambda = (Real*)malloc(lengthA*m*k*sizeof(Real));
		Real* rhoLambda = (Real*)malloc(m*m*k*sizeof(Real));
		Real* piLambda = (Real*)malloc(k*sizeof(Real));
		Real* LLF = (Real*)malloc((maxi+1)*sizeof(Real));
		Real* S = (Real*)malloc(lengthA*m*k*sizeof(Real));
		EMGLLF_core(phia,rhoInit,piInit,gamInit,mini,maxi,gamma,0.,Xa,Y,tau,
			phiLambda,rhoLambda,piLambda,LLF,S,
			n,lengthA,m,k);
		free(Xa);
		free(phia);
		free(LLF);
		free(S);

		//Assign results to current variables
		for (int j=0; j<lengthA; j++)
		{
			for (int mm=0; mm<m; mm++)
			{
				for (int r=0; r<k; r++)
					phi[ai4(a[j],mm,r,lambdaIndex,p,m,k,L)] = phiLambda[ai(j,mm,r,lengthA,m,k)];
			}
		}
		free(phiLambda);
		for (int u=0; u<m; u++)
		{
			for (int v=0; v<m; v++)
			{
				for (int r=0; r<k; r++)
					rho[ai4(u,v,r,lambdaIndex,m,m,k,L)] = rhoLambda[ai(u,v,r,m,m,k)];
			}
		}
		free(rhoLambda);
		for (int r=0; r<k; r++)
			pi[mi(r,lambdaIndex,k,L)] = piLambda[r];
		free(piLambda);

		int dimension = 0;
		int* b = (int*)malloc(m*sizeof(int));
		for (int j=0; j<p; j++)
		{
			//b = A2[j,2:dim(A2)[2],lambdaIndex] ; b = b[b!=0]
			int lengthB = 0;
			for (int mm=0; mm<m; mm++)
			{
				if (A2[ai(j,mm+1,lambdaIndex,p,m+1,L)] != 0)
					b[lengthB++] = A2[ai(j,mm+1,lambdaIndex,p,m+1,L)] - 1;
			}
			if (lengthB > 0)
			{
				//phi[A2[j,1,lambdaIndex],b,,lambdaIndex] = 0.
				for (int mm=0; mm<lengthB; mm++)
				{
					for (int r=0; r<k; r++)
						phi[ai4(A2[ai(j,0,lambdaIndex,p,m+1,L)]-1, b[mm], r, lambdaIndex, p, m, k, L)] = 0.;
				}
			}

			//c = A1[j,2:dim(A1)[2],lambdaIndex] ; dimension = dimension + sum(c!=0)
			for (int mm=0; mm<m; mm++)
			{
				if (A1[ai(j,mm+1,lambdaIndex,p,m+1,L)] != 0)
					dimension++;
			}
		}
		free(b);

		int signum;
		Real* densite = (Real*)calloc(L*n,sizeof(Real));
		Real sumLogDensit = 0.0;
		gsl_matrix* matrix = gsl_matrix_alloc(m, m);
		gsl_permutation* permutation = gsl_permutation_alloc(m);
		Real* YiRhoR = (Real*)malloc(m*sizeof(Real));
		Real* XiPhiR = (Real*)malloc(m*sizeof(Real));
		for (int i=0; i<n; i++)
		{
			for (int r=0; r<k; r++)
			{
				//compute det(rho(:,:,r,lambdaIndex)) [TODO: avoid re-computations]
				for (int u=0; u<m; u++)
				{
					for (int v=0; v<m; v++)
						matrix->data[u*m+v] = rho[ai4(u,v,r,lambdaIndex,m,m,k,L)];
				}
				gsl_linalg_LU_decomp(matrix, permutation, &signum);
				Real detRhoR = gsl_linalg_LU_det(matrix, signum);

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
				// NOTE: On peut remplacer X par Xa dans ce dernier calcul,
				// mais je ne sais pas si c'est intéressant ...

				// compute dotProduct < delta . delta >
				Real dotProduct = 0.0;
				for (int u=0; u<m; u++)
					dotProduct += (YiRhoR[u]-XiPhiR[u]) * (YiRhoR[u]-XiPhiR[u]);

				densite[mi(lambdaIndex,i,L,n)] +=
					(pi[mi(r,lambdaIndex,k,L)]*detRhoR/pow(sqrt(2.0*M_PI),m))*exp(-dotProduct/2.0);
			}
			sumLogDensit += log(densite[lambdaIndex*n+i]);
		}
		llh[mi(lambdaIndex,0,L,2)] = sumLogDensit;
		llh[mi(lambdaIndex,1,L,2)] = (dimension+m+1)*k-1;

		free(a);
		free(YiRhoR);
		free(XiPhiR);
		free(densite);
		gsl_matrix_free(matrix);
		gsl_permutation_free(permutation);
	}
	}
}
