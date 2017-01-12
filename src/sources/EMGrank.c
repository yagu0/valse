#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include "utils.h"

// Compute pseudo-inverse of a square matrix
static float* pinv(const float* matrix, int dim)
{
	gsl_matrix* U = gsl_matrix_alloc(dim,dim);
	gsl_matrix* V = gsl_matrix_alloc(dim,dim);
	gsl_vector* S = gsl_vector_alloc(dim);
	gsl_vector* work = gsl_vector_alloc(dim);
	float EPS = 1e-10; //threshold for singular value "== 0"
	
	//copy matrix into U
	copyArray(matrix, U->data, dim*dim);

	//U,S,V = SVD of matrix
	gsl_linalg_SV_decomp(U, V, S, work);
	gsl_vector_free(work);

	// Obtain pseudo-inverse by V*S^{-1}*t(U)
	float* inverse = (float*)malloc(dim*dim*sizeof(float));
	for (int i=0; i<dim; i++)
	{
		for (int ii=0; ii<dim; ii++)
		{
			float dotProduct = 0.0;
			for (int j=0; j<dim; j++)
				dotProduct += V->data[i*dim+j] * (S->data[j] > EPS ? 1.0/S->data[j] : 0.0) * U->data[ii*dim+j];
			inverse[i*dim+ii] = dotProduct;
		}
	}
	
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	return inverse;
}

// TODO: comment EMGrank purpose
void EMGrank_core(
	// IN parameters
	const float* Pi, // parametre de proportion
	const float* Rho, // parametre initial de variance renormalisé
	int mini, // nombre minimal d'itérations dans l'algorithme EM
	int maxi, // nombre maximal d'itérations dans l'algorithme EM
	const float* X, // régresseurs
	const float* Y, // réponse
	float tau, // seuil pour accepter la convergence
	const int* rank, // vecteur des rangs possibles
	// OUT parameters
	float* phi, // parametre de moyenne renormalisé, calculé par l'EM
	float* LLF, // log vraisemblance associé à cet échantillon, pour les valeurs estimées des paramètres
	// additional size parameters
	int n, // taille de l'echantillon
	int p, // nombre de covariables
	int m, // taille de Y (multivarié)
	int k) // nombre de composantes
{
	// Allocations, initializations
	float* Phi = (float*)calloc(p*m*k,sizeof(float));
	float* hatBetaR = (float*)malloc(p*m*sizeof(float));
	int signum;
	float invN = 1.0/n;
	int deltaPhiBufferSize = 20;
	float* deltaPhi = (float*)malloc(deltaPhiBufferSize*sizeof(float));
	int ite = 0;
	float sumDeltaPhi = 0.0;
	float* YiRhoR = (float*)malloc(m*sizeof(float));
	float* XiPhiR = (float*)malloc(m*sizeof(float));
	float* Xr = (float*)malloc(n*p*sizeof(float));
	float* Yr = (float*)malloc(n*m*sizeof(float));
	float* tXrXr = (float*)malloc(p*p*sizeof(float));
	float* tXrYr = (float*)malloc(p*m*sizeof(float));
	gsl_matrix* matrixM = gsl_matrix_alloc(p, m);
	gsl_matrix* matrixE = gsl_matrix_alloc(m, m);
	gsl_permutation* permutation = gsl_permutation_alloc(m);
	gsl_matrix* V = gsl_matrix_alloc(m,m);
	gsl_vector* S = gsl_vector_alloc(m);
	gsl_vector* work = gsl_vector_alloc(m);

	//Initialize class memberships (all elements in class 0; TODO: randomize ?)
	int* Z = (int*)calloc(n, sizeof(int));

	//Initialize phi to zero, because some M loops might exit before phi affectation
	zeroArray(phi, p*m*k);

	while (ite<mini || (ite<maxi && sumDeltaPhi>tau))
	{
		/////////////
		// Etape M //
		/////////////
		
		//M step: Mise à jour de Beta (et donc phi)
		for (int r=0; r<k; r++)
		{
			//Compute Xr = X(Z==r,:) and Yr = Y(Z==r,:)
			int cardClustR=0;
			for (int i=0; i<n; i++)
			{
				if (Z[i] == r)
				{
					for (int j=0; j<p; j++)
						Xr[mi(cardClustR,j,n,p)] = X[mi(i,j,n,p)];
					for (int j=0; j<m; j++)
						Yr[mi(cardClustR,j,n,m)] = Y[mi(i,j,n,m)];
					cardClustR++;
				}
			}
			if (cardClustR == 0)
				continue;

			//Compute tXrXr = t(Xr) * Xr
			for (int j=0; j<p; j++)
			{
				for (int jj=0; jj<p; jj++)
				{
					float dotProduct = 0.0;
					for (int u=0; u<cardClustR; u++)
						dotProduct += Xr[mi(u,j,n,p)] * Xr[mi(u,jj,n,p)];
					tXrXr[mi(j,jj,p,p)] = dotProduct;
				}
			}

			//Get pseudo inverse = (t(Xr)*Xr)^{-1}
			float* invTXrXr = pinv(tXrXr, p);
			
			// Compute tXrYr = t(Xr) * Yr
			for (int j=0; j<p; j++)
			{
				for (int jj=0; jj<m; jj++)
				{
					float dotProduct = 0.0;
					for (int u=0; u<cardClustR; u++)
						dotProduct += Xr[mi(u,j,n,p)] * Yr[mi(u,j,n,m)];
					tXrYr[mi(j,jj,p,m)] = dotProduct;
				}
			}

			//Fill matrixM with inverse * tXrYr = (t(Xr)*Xr)^{-1} * t(Xr) * Yr
			for (int j=0; j<p; j++)
			{
				for (int jj=0; jj<m; jj++)
				{
					float dotProduct = 0.0;
					for (int u=0; u<p; u++)
						dotProduct += invTXrXr[mi(j,u,p,p)] * tXrYr[mi(u,jj,p,m)];
					matrixM->data[j*m+jj] = dotProduct;
				}
			}
			free(invTXrXr);

			//U,S,V = SVD of (t(Xr)Xr)^{-1} * t(Xr) * Yr
			gsl_linalg_SV_decomp(matrixM, V, S, work);

			//Set m-rank(r) singular values to zero, and recompose
			//best rank(r) approximation of the initial product
			for (int j=rank[r]; j<m; j++)
				S->data[j] = 0.0;
			
			//[intermediate step] Compute hatBetaR = U * S * t(V)
			double* U = matrixM->data;
			for (int j=0; j<p; j++)
			{
				for (int jj=0; jj<m; jj++)
				{
					float dotProduct = 0.0;
					for (int u=0; u<m; u++)
						dotProduct += U[j*m+u] * S->data[u] * V->data[jj*m+u];
					hatBetaR[mi(j,jj,p,m)] = dotProduct;
				}
			}

			//Compute phi(:,:,r) = hatBetaR * Rho(:,:,r)
			for (int j=0; j<p; j++)
			{
				for (int jj=0; jj<m; jj++)
				{
					float dotProduct=0.0;
					for (int u=0; u<m; u++)
						dotProduct += hatBetaR[mi(j,u,p,m)] * Rho[ai(u,jj,r,m,m,k)];
					phi[ai(j,jj,r,p,m,k)] = dotProduct;
				}
		}
		}

		/////////////
		// Etape E //
		/////////////
		
		float sumLogLLF2 = 0.0;
		for (int i=0; i<n; i++)
		{
			float sumLLF1 = 0.0;
			float maxLogGamIR = -INFINITY;
			for (int r=0; r<k; r++)
			{
				//Compute
				//Gam(i,r) = Pi(r) * det(Rho(:,:,r)) * exp( -1/2 * (Y(i,:)*Rho(:,:,r) - X(i,:)...
				//*phi(:,:,r)) * transpose( Y(i,:)*Rho(:,:,r) - X(i,:)*phi(:,:,r) ) );
				//split in several sub-steps
				
				//compute det(Rho(:,:,r)) [TODO: avoid re-computations]
				for (int j=0; j<m; j++)
				{
					for (int jj=0; jj<m; jj++)
						matrixE->data[j*m+jj] = Rho[ai(j,jj,r,m,m,k)];
				}
				gsl_linalg_LU_decomp(matrixE, permutation, &signum);
				float detRhoR = gsl_linalg_LU_det(matrixE, signum);

				//compute Y(i,:)*Rho(:,:,r)
				for (int j=0; j<m; j++)
				{
					YiRhoR[j] = 0.0;
					for (int u=0; u<m; u++)
						YiRhoR[j] += Y[mi(i,u,n,m)] * Rho[ai(u,j,r,m,m,k)];
				}

				//compute X(i,:)*phi(:,:,r)
				for (int j=0; j<m; j++)
				{
					XiPhiR[j] = 0.0;
					for (int u=0; u<p; u++)
						XiPhiR[j] += X[mi(i,u,n,p)] * phi[ai(u,j,r,p,m,k)];
				}

				//compute dotProduct < Y(:,i)*rho(:,:,r)-X(i,:)*phi(:,:,r) . Y(:,i)*rho(:,:,r)-X(i,:)*phi(:,:,r) >
				float dotProduct = 0.0;
				for (int u=0; u<m; u++)
					dotProduct += (YiRhoR[u]-XiPhiR[u]) * (YiRhoR[u]-XiPhiR[u]);
				float logGamIR = log(Pi[r]) + log(detRhoR) - 0.5*dotProduct;

				//Z(i) = index of max (gam(i,:))
				if (logGamIR > maxLogGamIR)
				{
					Z[i] = r;
					maxLogGamIR = logGamIR;
				}
				sumLLF1 += exp(logGamIR) / pow(2*M_PI,m/2.0);
			}

			sumLogLLF2 += log(sumLLF1);
		}

		// Assign output variable LLF
		*LLF = -invN * sumLogLLF2;

		//newDeltaPhi = max(max((abs(phi-Phi))./(1+abs(phi))));
		float newDeltaPhi = 0.0;
		for (int j=0; j<p; j++)
		{
			for (int jj=0; jj<m; jj++)
			{
				for (int r=0; r<k; r++)
				{
					float tmpDist = fabs(phi[ai(j,jj,r,p,m,k)]-Phi[ai(j,jj,r,p,m,k)])
						/ (1.0+fabs(phi[ai(j,jj,r,p,m,k)]));
					if (tmpDist > newDeltaPhi)
						newDeltaPhi = tmpDist;
				}
			}
		}

		//update distance parameter to check algorithm convergence (delta(phi, Phi))
		//TODO: deltaPhi should be a linked list for perf.
		if (ite < deltaPhiBufferSize)
			deltaPhi[ite] = newDeltaPhi;
		else
		{
			sumDeltaPhi -= deltaPhi[0];
			for (int u=0; u<deltaPhiBufferSize-1; u++)
				deltaPhi[u] = deltaPhi[u+1];
			deltaPhi[deltaPhiBufferSize-1] = newDeltaPhi;
		}
		sumDeltaPhi += newDeltaPhi;

		// update other local variables
		for (int j=0; j<m; j++)
		{
			for (int jj=0; jj<p; jj++)
			{
				for (int r=0; r<k; r++)
					Phi[ai(j,jj,r,p,m,k)] = phi[ai(j,jj,r,p,m,k)];
			}
		}
		ite++;
	}

	//free memory
	free(hatBetaR);
	free(deltaPhi);
	free(Phi);
	gsl_matrix_free(matrixE);
	gsl_matrix_free(matrixM);
	gsl_permutation_free(permutation);
	gsl_vector_free(work);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	free(XiPhiR);
	free(YiRhoR);
	free(Xr);
	free(Yr);
	free(tXrXr);
	free(tXrYr);
	free(Z);
}
