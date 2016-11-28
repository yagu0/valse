#include "EMGLLF.h"
#include "constructionModelesLassoMLE.h"
#include <gsl/gsl_linalg.h>
#include <omp.h>
#include "omp_num_threads.h"

// TODO: comment on constructionModelesLassoMLE purpose
void constructionModelesLassoMLE(
	// IN parameters 
	const Real* phiInit, // parametre initial de moyenne renormalisé
	const Real* rhoInit, // parametre initial de variance renormalisé
	const Real* piInit,  // parametre initial des proportions
	const Real* gamInit, // paramètre initial des probabilités a posteriori de chaque échantillon
	Int mini,        // nombre minimal d'itérations dans l'algorithme EM	
	Int maxi,        // nombre maximal d'itérations dans l'algorithme EM
	Real gamma,    // valeur de gamma : puissance des proportions dans la pénalisation pour un Lasso adaptatif
	const Real* glambda, // valeur des paramètres de régularisation du Lasso
	const Real* X,       // régresseurs
	const Real* Y,       // réponse
	Real seuil,    // seuil pour prendre en compte une variable
	Real tau,      // seuil pour accepter la convergence
	const Int* A1,         // matrice des coefficients des parametres selectionnes
	const Int* A2,         // matrice des coefficients des parametres non selectionnes
	// OUT parameters
	Real* phi,            // estimateur ainsi calculé par le Lasso
    Real* rho,            // estimateur ainsi calculé par le Lasso
	Real* pi,             // estimateur ainsi calculé par le Lasso
    Real* lvraisemblance, // estimateur ainsi calculé par le Lasso
	// additional size parameters
	mwSize n,                 // taille de l'echantillon                
	mwSize p,                 // nombre de covariables
	mwSize m,                 // taille de Y (multivarié)
	mwSize k,                 // nombre de composantes
	mwSize L)                 // taille de glambda
{
	//preparation: phi = 0
	for (mwSize u=0; u<p*m*k*L; u++)
		phi[u] = 0.0;
	
	//initiate parallel section
	mwSize lambdaIndex;
	omp_set_num_threads(OMP_NUM_THREADS);
	#pragma omp parallel default(shared) private(lambdaIndex)
	{
	#pragma omp for schedule(dynamic,CHUNK_SIZE) nowait
	for (lambdaIndex=0; lambdaIndex<L; lambdaIndex++)
	{
		//~ a = A1(:,1,lambdaIndex);
		//~ a(a==0) = [];
		Int* a = (Int*)malloc(p*sizeof(Int));
		mwSize lengthA = 0;
		for (mwSize j=0; j<p; j++)
		{
			if (A1[j*(m+1)*L+0*L+lambdaIndex] != 0)
				a[lengthA++] = A1[j*(m+1)*L+0*L+lambdaIndex] - 1;
		}
		if (lengthA == 0)
			continue;
		
		//Xa = X(:,a)
		Real* Xa = (Real*)malloc(n*lengthA*sizeof(Real));
		for (mwSize i=0; i<n; i++)
		{
			for (mwSize j=0; j<lengthA; j++)
				Xa[i*lengthA+j] = X[i*p+a[j]];
		}
		
		//phia = phiInit(a,:,:)
		Real* phia = (Real*)malloc(lengthA*m*k*sizeof(Real));
		for (mwSize j=0; j<lengthA; j++)
		{
			for (mwSize mm=0; mm<m; mm++)
			{
				for (mwSize r=0; r<k; r++)
					phia[j*m*k+mm*k+r] = phiInit[a[j]*m*k+mm*k+r];
			}
		}
		
		//[phiLambda,rhoLambda,piLambda,~,~] = EMGLLF(...
		//	phiInit(a,:,:),rhoInit,piInit,gamInit,mini,maxi,gamma,0,X(:,a),Y,tau);
		Real* phiLambda = (Real*)malloc(lengthA*m*k*sizeof(Real));
		Real* rhoLambda = (Real*)malloc(m*m*k*sizeof(Real));
		Real* piLambda = (Real*)malloc(k*sizeof(Real));
		Real* LLF = (Real*)malloc((maxi+1)*sizeof(Real));
		Real* S = (Real*)malloc(lengthA*m*k*sizeof(Real));
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
		for (mwSize j=0; j<lengthA; j++)
		{
			for (mwSize mm=0; mm<m; mm++)
			{
				for (mwSize r=0; r<k; r++)
					phi[a[j]*m*k*L+mm*k*L+r*L+lambdaIndex] = phiLambda[j*m*k+mm*k+r];
			}
		}
		free(phiLambda);
		//~ rho(:,:,:,lambdaIndex) = rhoLambda;
		for (mwSize u=0; u<m; u++)
		{
			for (mwSize v=0; v<m; v++)
			{
				for (mwSize r=0; r<k; r++)
					rho[u*m*k*L+v*k*L+r*L+lambdaIndex] = rhoLambda[u*m*k+v*k+r];
			}
		}
		free(rhoLambda);
		//~ pi(:,lambdaIndex) = piLambda;
		for (mwSize r=0; r<k; r++)
			pi[r*L+lambdaIndex] = piLambda[r];
		free(piLambda);
		
		mwSize dimension = 0;
		Int* b = (Int*)malloc(m*sizeof(Int));
		for (mwSize j=0; j<p; j++)
		{
			//~ b = A2(j,2:end,lambdaIndex);
			//~ b(b==0) = [];
			mwSize lengthB = 0;
			for (mwSize mm=0; mm<m; mm++)
			{
				if (A2[j*(m+1)*L+(mm+1)*L+lambdaIndex] != 0)
					b[lengthB++] = A2[j*(m+1)*L+(mm+1)*L+lambdaIndex] - 1;
			}
			//~ if length(b) > 0
				//~ phi(A2(j,1,lambdaIndex),b,:,lambdaIndex) = 0.0;
			//~ end
			if (lengthB > 0)
			{
				for (mwSize mm=0; mm<lengthB; mm++)
				{
					for (mwSize r=0; r<k; r++)
						phi[(A2[j*(m+1)*L+0*L+lambdaIndex]-1)*m*k*L + b[mm]*k*L + r*L + lambdaIndex] = 0.0;
				}
			}
			
			//~ c = A1(j,2:end,lambdaIndex);
			//~ c(c==0) = [];
			//~ dimension = dimension + length(c);
			for (mwSize mm=0; mm<m; mm++)
			{
				if (A1[j*(m+1)*L+(mm+1)*L+lambdaIndex] != 0)
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
		for (mwSize i=0; i<n; i++)
		{
			//~ for r=1:k
				//~ delta = Y(i,:)*rho(:,:,r,lambdaIndex) - (X(i,a)*(phi(a,:,r,lambdaIndex)));
				//~ densite(i,lambdaIndex) = densite(i,lambdaIndex) +...
					//~ pi(r,lambdaIndex)*det(rho(:,:,r,lambdaIndex))/(sqrt(2*PI))^m*exp(-dot(delta,delta)/2.0);
			//~ end
			for (mwSize r=0; r<k; r++)
			{
				//compute det(rho(:,:,r,lambdaIndex)) [TODO: avoid re-computations]
				for (mwSize u=0; u<m; u++)
				{
					for (mwSize v=0; v<m; v++)
						matrix->data[u*m+v] = rho[u*m*k*L+v*k*L+r*L+lambdaIndex];
				}
				gsl_linalg_LU_decomp(matrix, permutation, &signum);
				Real detRhoR = gsl_linalg_LU_det(matrix, signum);
				
				//compute Y(i,:)*rho(:,:,r,lambdaIndex)
				for (mwSize u=0; u<m; u++)
				{
					YiRhoR[u] = 0.0;
					for (mwSize v=0; v<m; v++)
						YiRhoR[u] += Y[i*m+v] * rho[v*m*k*L+u*k*L+r*L+lambdaIndex];
				}
				
				//compute X(i,a)*phi(a,:,r,lambdaIndex)
				for (mwSize u=0; u<m; u++)
				{
					XiPhiR[u] = 0.0;
					for (mwSize v=0; v<lengthA; v++)
						XiPhiR[u] += X[i*p+a[v]] * phi[a[v]*m*k*L+u*k*L+r*L+lambdaIndex];
				}
                // On peut remplacer X par Xa dans ce dernier calcul, mais je ne sais pas si c'est intéressant ...
				
				// compute dotProduct < delta . delta >
				Real dotProduct = 0.0;
				for (mwSize u=0; u<m; u++)
					dotProduct += (YiRhoR[u]-XiPhiR[u]) * (YiRhoR[u]-XiPhiR[u]);
				
				densite[lambdaIndex*n+i] += (pi[r*L+lambdaIndex]*detRhoR/pow(sqrt(2.0*M_PI),m))*exp(-dotProduct/2.0);
			}			
			sumLogDensit += log(densite[lambdaIndex*n+i]);
		}
		lvraisemblance[lambdaIndex*2+0] = sumLogDensit;
		lvraisemblance[lambdaIndex*2+1] = (dimension+m+1)*k-1;
	
		free(a);
		free(YiRhoR);
		free(XiPhiR);
		free(densite);
		gsl_matrix_free(matrix);
		gsl_permutation_free(permutation);
	}
	}
}
