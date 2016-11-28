#include "EMGrank.h"
#include <gsl/gsl_linalg.h>

// Compute pseudo-inverse of a square matrix
static Real* pinv(const Real* matrix, mwSize dim)
{
	gsl_matrix* U = gsl_matrix_alloc(dim,dim);
	gsl_matrix* V = gsl_matrix_alloc(dim,dim);
	gsl_vector* S = gsl_vector_alloc(dim);
	gsl_vector* work = gsl_vector_alloc(dim);
	Real EPS = 1e-10; //threshold for singular value "== 0"
	
	//copy matrix into U
	for (mwSize i=0; i<dim*dim; i++)
		U->data[i] = matrix[i];
	
	//U,S,V = SVD of matrix
	gsl_linalg_SV_decomp(U, V, S, work);
	gsl_vector_free(work);
	
	// Obtain pseudo-inverse by V*S^{-1}*t(U)
	Real* inverse = (Real*)malloc(dim*dim*sizeof(Real));
	for (mwSize i=0; i<dim; i++)
	{
		for (mwSize ii=0; ii<dim; ii++)
		{
			Real dotProduct = 0.0;
			for (mwSize j=0; j<dim; j++)
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
void EMGrank(
	// IN parameters
	const Real* Pi,          // parametre de proportion
	const Real* Rho,         // parametre initial de variance renormalisé
	Int mini,            // nombre minimal d'itérations dans l'algorithme EM        
	Int maxi,            // nombre maximal d'itérations dans l'algorithme EM
	const Real* X,           // régresseurs
	const Real* Y,           // réponse
	Real tau,          // seuil pour accepter la convergence
	const Int* rank,           // vecteur des rangs possibles
	// OUT parameters
	Real* phi,         // parametre de moyenne renormalisé, calculé par l'EM
	Real* LLF,         // log vraisemblance associé à cet échantillon, pour les valeurs estimées des paramètres
	// additional size parameters
	mwSize n,               // taille de l'echantillon                
	mwSize p,               // nombre de covariables
	mwSize m,               // taille de Y (multivarié)
	mwSize k)              // nombre de composantes
{
	// Allocations, initializations
	Real* Phi = (Real*)calloc(p*m*k,sizeof(Real));
	Real* hatBetaR = (Real*)malloc(p*m*sizeof(Real));
	int signum;
	Real invN = 1.0/n;
	int deltaPhiBufferSize = 20;
	Real* deltaPhi = (Real*)malloc(deltaPhiBufferSize*sizeof(Real));
	mwSize ite = 0;
	Real sumDeltaPhi = 0.0;
	Real* YiRhoR = (Real*)malloc(m*sizeof(Real));
	Real* XiPhiR = (Real*)malloc(m*sizeof(Real));
	Real* Xr = (Real*)malloc(n*p*sizeof(Real));
	Real* Yr = (Real*)malloc(n*m*sizeof(Real));
	Real* tXrXr = (Real*)malloc(p*p*sizeof(Real));
	Real* tXrYr = (Real*)malloc(p*m*sizeof(Real));
	gsl_matrix* matrixM = gsl_matrix_alloc(p, m);        
	gsl_matrix* matrixE = gsl_matrix_alloc(m, m);
	gsl_permutation* permutation = gsl_permutation_alloc(m);
	gsl_matrix* V = gsl_matrix_alloc(m,m);
	gsl_vector* S = gsl_vector_alloc(m);
	gsl_vector* work = gsl_vector_alloc(m);

	//Initialize class memberships (all elements in class 0; TODO: randomize ?)
	Int* Z = (Int*)calloc(n, sizeof(Int));

	//Initialize phi to zero, because some M loops might exit before phi affectation
	for (mwSize i=0; i<p*m*k; i++)
		phi[i] = 0.0;

	while (ite<mini || (ite<maxi && sumDeltaPhi>tau))
	{		
		/////////////
		// Etape M //
		/////////////
		
		//M step: Mise à jour de Beta (et donc phi)
		for (mwSize r=0; r<k; r++)
		{
			//Compute Xr = X(Z==r,:) and Yr = Y(Z==r,:)
			mwSize cardClustR=0;
			for (mwSize i=0; i<n; i++)
			{
				if (Z[i] == r)
				{
					for (mwSize j=0; j<p; j++)
						Xr[cardClustR*p+j] = X[i*p+j];
					for (mwSize j=0; j<m; j++)
						Yr[cardClustR*m+j] = Y[i*m+j];
					cardClustR++;
				}
			}
			if (cardClustR == 0) 
				continue;

			//Compute tXrXr = t(Xr) * Xr
			for (mwSize j=0; j<p; j++)
			{
				for (mwSize jj=0; jj<p; jj++)
				{
					Real dotProduct = 0.0;
					for (mwSize u=0; u<cardClustR; u++)
						dotProduct += Xr[u*p+j] * Xr[u*p+jj];
					tXrXr[j*p+jj] = dotProduct;
				}
			}

			//Get pseudo inverse = (t(Xr)*Xr)^{-1}
			Real* invTXrXr = pinv(tXrXr, p);
			
			// Compute tXrYr = t(Xr) * Yr
			for (mwSize j=0; j<p; j++)
			{
				for (mwSize jj=0; jj<m; jj++)
				{
					Real dotProduct = 0.0;
					for (mwSize u=0; u<cardClustR; u++)
						dotProduct += Xr[u*p+j] * Yr[u*m+jj];
					tXrYr[j*m+jj] = dotProduct;
				}
			}

			//Fill matrixM with inverse * tXrYr = (t(Xr)*Xr)^{-1} * t(Xr) * Yr
			for (mwSize j=0; j<p; j++)
			{
				for (mwSize jj=0; jj<m; jj++)
				{
					Real dotProduct = 0.0;
					for (mwSize u=0; u<p; u++)
						dotProduct += invTXrXr[j*p+u] * tXrYr[u*m+jj];
					matrixM->data[j*m+jj] = dotProduct;	
				}
			}
			free(invTXrXr);
			
			//U,S,V = SVD of (t(Xr)Xr)^{-1} * t(Xr) * Yr
			gsl_linalg_SV_decomp(matrixM, V, S, work);
			
			//Set m-rank(r) singular values to zero, and recompose 
			//best rank(r) approximation of the initial product
			for (mwSize j=rank[r]; j<m; j++)
				S->data[j] = 0.0;
			
			//[intermediate step] Compute hatBetaR = U * S * t(V)
			Real* U = matrixM->data;
			for (mwSize j=0; j<p; j++)
			{
				for (mwSize jj=0; jj<m; jj++)
				{
					Real dotProduct = 0.0;
					for (mwSize u=0; u<m; u++) 
						dotProduct += U[j*m+u] * S->data[u] * V->data[jj*m+u];
					hatBetaR[j*m+jj] = dotProduct;
				}
			}
 
			//Compute phi(:,:,r) = hatBetaR * Rho(:,:,r)
			for (mwSize j=0; j<p; j++)
			{
				for (mwSize jj=0; jj<m; jj++)
				{
					Real dotProduct=0.0;
					for (mwSize u=0; u<m; u++)
						dotProduct += hatBetaR[j*m+u] * Rho[u*m*k+jj*k+r];
					phi[j*m*k+jj*k+r] = dotProduct;
				}
		    }
		}

		/////////////
		// Etape E //
		/////////////
		
		Real sumLogLLF2 = 0.0;
		for (mwSize i=0; i<n; i++)
		{
			Real sumLLF1 = 0.0;
			Real maxLogGamIR = -INFINITY;
			for (mwSize r=0; r<k; r++)
			{
				//Compute
				//Gam(i,r) = Pi(r) * det(Rho(:,:,r)) * exp( -1/2 * (Y(i,:)*Rho(:,:,r) - X(i,:)...
				//    *phi(:,:,r)) * transpose( Y(i,:)*Rho(:,:,r) - X(i,:)*phi(:,:,r) ) );
				//split in several sub-steps
				
				//compute det(Rho(:,:,r)) [TODO: avoid re-computations]
				for (mwSize j=0; j<m; j++)
				{
					for (mwSize jj=0; jj<m; jj++)
						matrixE->data[j*m+jj] = Rho[j*m*k+jj*k+r];
				}
				gsl_linalg_LU_decomp(matrixE, permutation, &signum);
				Real detRhoR = gsl_linalg_LU_det(matrixE, signum);
				
				//compute Y(i,:)*Rho(:,:,r)
				for (mwSize j=0; j<m; j++)
				{
					YiRhoR[j] = 0.0;
					for (mwSize u=0; u<m; u++)
						YiRhoR[j] += Y[i*m+u] * Rho[u*m*k+j*k+r];
				}
				
				//compute X(i,:)*phi(:,:,r)
				for (mwSize j=0; j<m; j++)
				{
					XiPhiR[j] = 0.0;
					for (mwSize u=0; u<p; u++)
						XiPhiR[j] += X[i*p+u] * phi[u*m*k+j*k+r];
				}
				
				//compute dotProduct < Y(:,i)*rho(:,:,r)-X(i,:)*phi(:,:,r) . Y(:,i)*rho(:,:,r)-X(i,:)*phi(:,:,r) >
				Real dotProduct = 0.0;
				for (mwSize u=0; u<m; u++)
					dotProduct += (YiRhoR[u]-XiPhiR[u]) * (YiRhoR[u]-XiPhiR[u]);
				Real logGamIR = log(Pi[r]) + log(detRhoR) - 0.5*dotProduct;
				
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
		Real newDeltaPhi = 0.0;
		for (mwSize j=0; j<p; j++)
		{
			for (mwSize jj=0; jj<m; jj++)
			{
				for (mwSize r=0; r<k; r++)
				{
					Real tmpDist = fabs(phi[j*m*k+jj*k+r]-Phi[j*m*k+jj*k+r]) 
						/ (1.0+fabs(phi[j*m*k+jj*k+r]));
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
		for (mwSize j=0; j<m; j++)
		{
			for (mwSize jj=0; jj<p; jj++)
			{
				for (mwSize r=0; r<k; r++)
					Phi[j*m*k+jj*k+r] = phi[j*m*k+jj*k+r];
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
