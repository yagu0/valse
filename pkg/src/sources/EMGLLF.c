#include "utils.h"
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

// TODO: don't recompute indexes ai(...) and mi(...) when possible
void EMGLLF_core(
	// IN parameters
	const Real* phiInit, // parametre initial de moyenne renormalisé
	const Real* rhoInit, // parametre initial de variance renormalisé
	const Real* piInit,	 // parametre initial des proportions
	const Real* gamInit, // paramètre initial des probabilités a posteriori de chaque échantillon
	int mini, // nombre minimal d'itérations dans l'algorithme EM
	int maxi, // nombre maximal d'itérations dans l'algorithme EM
	Real gamma, // puissance des proportions dans la pénalisation pour un Lasso adaptatif
	Real lambda, // valeur du paramètre de régularisation du Lasso
	const Real* X, // régresseurs
	const Real* Y, // réponse
	Real tau, // seuil pour accepter la convergence
	// OUT parameters (all pointers, to be modified)
	Real* phi, // parametre de moyenne renormalisé, calculé par l'EM
	Real* rho, // parametre de variance renormalisé, calculé par l'EM
	Real* pi, // parametre des proportions renormalisé, calculé par l'EM
	Real* LLF, // log vraisemblance associée à cet échantillon, pour les valeurs estimées des paramètres
	Real* S,
	// additional size parameters
	int n, // nombre d'echantillons
	int p, // nombre de covariables
	int m, // taille de Y (multivarié)
	int k) // nombre de composantes dans le mélange
{
	//Initialize outputs
	copyArray(phiInit, phi, p*m*k);
	copyArray(rhoInit, rho, m*m*k);
	copyArray(piInit, pi, k);
	zeroArray(LLF, maxi);
	//S is already allocated, and doesn't need to be 'zeroed'

	//Other local variables: same as in R
	Real* gam = (Real*)malloc(n*k*sizeof(Real));
	copyArray(gamInit, gam, n*k);
	Real* Gram2 = (Real*)malloc(p*p*k*sizeof(Real));
	Real* ps2 = (Real*)malloc(p*m*k*sizeof(Real));
	Real* b = (Real*)malloc(k*sizeof(Real));
	Real* X2 = (Real*)malloc(n*p*k*sizeof(Real));
	Real* Y2 = (Real*)malloc(n*m*k*sizeof(Real));
	Real dist = 0.;
	Real dist2 = 0.;
	int ite = 0;
	Real* pi2 = (Real*)malloc(k*sizeof(Real));
	Real* ps = (Real*)malloc(m*k*sizeof(Real));
	Real* nY2 = (Real*)malloc(m*k*sizeof(Real));
	Real* ps1 = (Real*)malloc(n*m*k*sizeof(Real));
	Real* Gam = (Real*)malloc(n*k*sizeof(Real));
	const Real EPS = 1e-15;
	// Additional (not at this place, in R file)
	Real* gam2 = (Real*)malloc(k*sizeof(Real));
	Real* nY21 = (Real*)malloc(n*m*k*sizeof(Real));
	Real* sqNorm2 = (Real*)malloc(k*sizeof(Real));
	gsl_matrix* matrix = gsl_matrix_alloc(m, m);
	gsl_permutation* permutation = gsl_permutation_alloc(m);
	Real* YiRhoR = (Real*)malloc(m*sizeof(Real));
	Real* XiPhiR = (Real*)malloc(m*sizeof(Real));
	const Real gaussConstM = pow(2.*M_PI,m/2.);
	Real* Phi = (Real*)malloc(p*m*k*sizeof(Real));
	Real* Rho = (Real*)malloc(m*m*k*sizeof(Real));
	Real* Pi = (Real*)malloc(k*sizeof(Real));

	while (ite < mini || (ite < maxi && (dist >= tau || dist2 >= sqrt(tau))))
	{
		copyArray(phi, Phi, p*m*k);
		copyArray(rho, Rho, m*m*k);
		copyArray(pi, Pi, k);

		// Calculs associés a Y et X
		for (int r=0; r<k; r++)
		{
			for (int mm=0; mm<m; mm++)
			{
				//Y2[,mm,r] = sqrt(gam[,r]) * Y[,mm]
				for (int u=0; u<n; u++)
					Y2[ai(u,mm,r,n,m,k)] = sqrt(gam[mi(u,r,n,k)]) * Y[mi(u,mm,m,n)];
			}
			for (int i=0; i<n; i++)
			{
				//X2[i,,r] = sqrt(gam[i,r]) * X[i,]
				for (int u=0; u<p; u++)
					X2[ai(i,u,r,n,p,k)] = sqrt(gam[mi(i,r,n,k)]) * X[mi(i,u,n,p)];
			}
			for (int mm=0; mm<m; mm++)
			{
				//ps2[,mm,r] = crossprod(X2[,,r],Y2[,mm,r])
				for (int u=0; u<p; u++)
				{
					Real dotProduct = 0.;
					for (int v=0; v<n; v++)
						dotProduct += X2[ai(v,u,r,n,p,k)] * Y2[ai(v,mm,r,n,m,k)];
					ps2[ai(u,mm,r,p,m,k)] = dotProduct;
				}
			}
			for (int j=0; j<p; j++)
			{
				for (int s=0; s<p; s++)
				{
					//Gram2[j,s,r] = crossprod(X2[,j,r], X2[,s,r])
					Real dotProduct = 0.;
					for (int u=0; u<n; u++)
						dotProduct += X2[ai(u,j,r,n,p,k)] * X2[ai(u,s,r,n,p,k)];
					Gram2[ai(j,s,r,p,p,k)] = dotProduct;
				}
			}
		}

		/////////////
		// Etape M //
		/////////////

		// Pour pi
		for (int r=0; r<k; r++)
		{
			//b[r] = sum(abs(phi[,,r]))
			Real sumAbsPhi = 0.;
			for (int u=0; u<p; u++)
				for (int v=0; v<m; v++)
					sumAbsPhi += fabs(phi[ai(u,v,r,p,m,k)]);
			b[r] = sumAbsPhi;
		}
		//gam2 = colSums(gam)
		for (int u=0; u<k; u++)
		{
			Real sumOnColumn = 0.;
			for (int v=0; v<n; v++)
				sumOnColumn += gam[mi(v,u,n,k)];
			gam2[u] = sumOnColumn;
		}
		//a = sum(gam %*% log(pi))
		Real a = 0.;
		for (int u=0; u<n; u++)
		{
			Real dotProduct = 0.;
			for (int v=0; v<k; v++)
				dotProduct += gam[mi(u,v,n,k)] * log(pi[v]);
			a += dotProduct;
		}

		//tant que les proportions sont negatives
		int kk = 0;
		int pi2AllPositive = 0;
		Real invN = 1./n;
		while (!pi2AllPositive)
		{
			//pi2 = pi + 0.1^kk * ((1/n)*gam2 - pi)
			Real pow_01_kk = pow(0.1,kk);
			for (int r=0; r<k; r++)
				pi2[r] = pi[r] + pow_01_kk * (invN*gam2[r] - pi[r]);
			//pi2AllPositive = all(pi2 >= 0)
			pi2AllPositive = 1;
			for (int r=0; r<k; r++)
			{
				if (pi2[r] < 0)
				{
					pi2AllPositive = 0;
					break;
				}
			}
			kk++;
		}

		//(pi.^gamma)*b
		Real piPowGammaDotB = 0.;
		for (int v=0; v<k; v++)
			piPowGammaDotB += pow(pi[v],gamma) * b[v];
		//(pi2.^gamma)*b
		Real pi2PowGammaDotB = 0.;
		for (int v=0; v<k; v++)
			pi2PowGammaDotB += pow(pi2[v],gamma) * b[v];
		//transpose(gam2)*log(pi2)
		Real prodGam2logPi2 = 0.;
		for (int v=0; v<k; v++)
			prodGam2logPi2 += gam2[v] * log(pi2[v]);
		//t(m) la plus grande valeur dans la grille O.1^k tel que ce soit décroissante ou constante
		while (-invN*a + lambda*piPowGammaDotB < -invN*prodGam2logPi2 + lambda*pi2PowGammaDotB
			&& kk<1000)
		{
			Real pow_01_kk = pow(0.1,kk);
			//pi2 = pi + 0.1^kk * (1/n*gam2 - pi)
			for (int v=0; v<k; v++)
				pi2[v] = pi[v] + pow_01_kk * (invN*gam2[v] - pi[v]);
			//pi2 was updated, so we recompute pi2PowGammaDotB and prodGam2logPi2
			pi2PowGammaDotB = 0.;
			for (int v=0; v<k; v++)
				pi2PowGammaDotB += pow(pi2[v],gamma) * b[v];
			prodGam2logPi2 = 0.;
			for (int v=0; v<k; v++)
				prodGam2logPi2 += gam2[v] * log(pi2[v]);
			kk++;
		}
		Real t = pow(0.1,kk);
		//sum(pi + t*(pi2-pi))
		Real sumPiPlusTbyDiff = 0.;
		for (int v=0; v<k; v++)
			sumPiPlusTbyDiff += (pi[v] + t*(pi2[v] - pi[v]));
		//pi = (pi + t*(pi2-pi)) / sum(pi + t*(pi2-pi))
		for (int v=0; v<k; v++)
			pi[v] = (pi[v] + t*(pi2[v] - pi[v])) / sumPiPlusTbyDiff;

		//Pour phi et rho
		for (int r=0; r<k; r++)
		{
			for (int mm=0; mm<m; mm++)
			{
				for (int i=0; i<n; i++)
				{
					//< X2(i,:,r) , phi(:,mm,r) >
					Real dotProduct = 0.;
					for (int u=0; u<p; u++)
						dotProduct += X2[ai(i,u,r,n,p,k)] * phi[ai(u,mm,r,p,m,k)];
					//ps1[i,mm,r] = Y2[i,mm,r] * sum(X2[i,,r] * phi[,mm,r])
					ps1[ai(i,mm,r,n,m,k)] = Y2[ai(i,mm,r,n,m,k)] * dotProduct;
					nY21[ai(i,mm,r,n,m,k)] = Y2[ai(i,mm,r,n,m,k)] * Y2[ai(i,mm,r,n,m,k)];
				}
				//ps[mm,r] = sum(ps1[,mm,r])
				Real sumPs1 = 0.;
				for (int u=0; u<n; u++)
					sumPs1 += ps1[ai(u,mm,r,n,m,k)];
				ps[mi(mm,r,m,k)] = sumPs1;
				//nY2[mm,r] = sum(nY21[,mm,r])
				Real sumNy21 = 0.;
				for (int u=0; u<n; u++)
					sumNy21 += nY21[ai(u,mm,r,n,m,k)];
				nY2[mi(mm,r,m,k)] = sumNy21;
				//rho[mm,mm,r] = (ps[mm,r]+sqrt(ps[mm,r]^2+4*nY2[mm,r]*(gam2[r]))) / (2*nY2[mm,r])
				rho[ai(mm,mm,r,m,m,k)] = ( ps[mi(mm,r,m,k)] + sqrt( ps[mi(mm,r,m,k)]*ps[mi(mm,r,m,k)]
					+ 4*nY2[mi(mm,r,m,k)] * gam2[r] ) ) / (2*nY2[mi(mm,r,m,k)]);
			}
		}
		for (int r=0; r<k; r++)
		{
			for (int j=0; j<p; j++)
			{
				for (int mm=0; mm<m; mm++)
				{
					//sum(phi[-j,mm,r] * Gram2[j, setdiff(1:p,j),r])
					Real dotPhiGram2 = 0.0;
					for (int u=0; u<p; u++)
					{
						if (u != j)
							dotPhiGram2 += phi[ai(u,mm,r,p,m,k)] * Gram2[ai(j,u,r,p,p,k)];
					}
					//S[j,mm,r] = -rho[mm,mm,r]*ps2[j,mm,r] + sum(phi[-j,mm,r] * Gram2[j, setdiff(1:p,j),r])
					S[ai(j,mm,r,p,m,k)] = -rho[ai(mm,mm,r,m,m,k)] * ps2[ai(j,mm,r,p,m,k)] + dotPhiGram2;
					Real pow_pir_gamma = pow(pi[r],gamma);
					if (fabs(S[ai(j,mm,r,p,m,k)]) <= n*lambda*pow_pir_gamma)
						phi[ai(j,mm,r,p,m,k)] = 0;
					else if (S[ai(j,mm,r,p,m,k)] > n*lambda*pow_pir_gamma)
					{
						phi[ai(j,mm,r,p,m,k)] = (n*lambda*pow_pir_gamma - S[ai(j,mm,r,p,m,k)])
							/ Gram2[ai(j,j,r,p,p,k)];
					}
					else
					{
						phi[ai(j,mm,r,p,m,k)] = -(n*lambda*pow_pir_gamma + S[ai(j,mm,r,p,m,k)])
							/ Gram2[ai(j,j,r,p,p,k)];
					}
				}
			}
		}

		/////////////
		// Etape E //
		/////////////

		int signum;
		Real sumLogLLF2 = 0.;
		for (int i=0; i<n; i++)
		{
			for (int r=0; r<k; r++)
			{
				//compute Y[i,]%*%rho[,,r]
				for (int u=0; u<m; u++)
				{
					YiRhoR[u] = 0.;
					for (int v=0; v<m; v++)
						YiRhoR[u] += Y[mi(i,v,n,m)] * rho[ai(v,u,r,m,m,k)];
				}

				//compute X(i,:)*phi(:,:,r)
				for (int u=0; u<m; u++)
				{
					XiPhiR[u] = 0.;
					for (int v=0; v<p; v++)
						XiPhiR[u] += X[mi(i,v,n,p)] * phi[ai(v,u,r,p,m,k)];
				}

				//compute sq norm || Y(:,i)*rho(:,:,r)-X(i,:)*phi(:,:,r) ||_2^2
				sqNorm2[r] = 0.;
				for (int u=0; u<m; u++)
					sqNorm2[r] += (YiRhoR[u]-XiPhiR[u]) * (YiRhoR[u]-XiPhiR[u]);
			}

			Real sumLLF1 = 0.;
			Real sumGamI = 0.;
			for (int r=0; r<k; r++)
			{
				//compute det(rho[,,r]) [TODO: avoid re-computations]
				for (int u=0; u<m; u++)
				{
					for (int v=0; v<m; v++)
						matrix->data[u*m+v] = rho[ai(u,v,r,m,m,k)];
				}
				gsl_linalg_LU_decomp(matrix, permutation, &signum);
				Real detRhoR = gsl_linalg_LU_det(matrix, signum);
				Gam[mi(i,r,n,k)] = pi[r] * exp(-.5*sqNorm2[r]) * detRhoR;
				sumLLF1 += Gam[mi(i,r,n,k)] / gaussConstM;
				sumGamI += Gam[mi(i,r,n,k)];
			}
			sumLogLLF2 += log(sumLLF1);
			for (int r=0; r<k; r++)
			{
				//gam[i,] = Gam[i,] / sumGamI
				gam[mi(i,r,n,k)] = sumGamI > EPS ? Gam[mi(i,r,n,k)] / sumGamI : 0.;
			}
		}

		//sumPen = sum(pi^gamma * b)
		Real sumPen = 0.;
		for (int r=0; r<k; r++)
			sumPen += pow(pi[r],gamma) * b[r];
		//LLF[ite] = -sumLogLLF2/n + lambda*sumPen
		LLF[ite] = -invN * sumLogLLF2 + lambda * sumPen;
		dist = ite==0 ? LLF[ite] : (LLF[ite] - LLF[ite-1]) / (1. + fabs(LLF[ite]));

		//Dist1 = max( abs(phi-Phi) / (1+abs(phi)) )
		Real Dist1 = 0.;
		for (int u=0; u<p; u++)
		{
			for (int v=0; v<m; v++)
			{
				for (int w=0; w<k; w++)
				{
					Real tmpDist = fabs(phi[ai(u,v,w,p,m,k)]-Phi[ai(u,v,w,p,m,k)])
						/ (1.+fabs(phi[ai(u,v,w,p,m,k)]));
					if (tmpDist > Dist1)
						Dist1 = tmpDist;
				}
			}
		}
		//Dist2 = max( (abs(rho-Rho)) / (1+abs(rho)) )
		Real Dist2 = 0.;
		for (int u=0; u<m; u++)
		{
			for (int v=0; v<m; v++)
			{
				for (int w=0; w<k; w++)
				{
					Real tmpDist = fabs(rho[ai(u,v,w,m,m,k)]-Rho[ai(u,v,w,m,m,k)])
						/ (1.+fabs(rho[ai(u,v,w,m,m,k)]));
					if (tmpDist > Dist2)
						Dist2 = tmpDist;
				}
			}
		}
		//Dist3 = max( (abs(pi-Pi)) / (1+abs(Pi)))
		Real Dist3 = 0.;
		for (int u=0; u<n; u++)
		{
			for (int v=0; v<k; v++)
			{
				Real tmpDist = fabs(pi[v]-Pi[v]) / (1.+fabs(pi[v]));
				if (tmpDist > Dist3)
					Dist3 = tmpDist;
			}
		}
		//dist2=max([max(Dist1),max(Dist2),max(Dist3)]);
		dist2 = Dist1;
		if (Dist2 > dist2)
			dist2 = Dist2;
		if (Dist3 > dist2)
			dist2 = Dist3;

		ite++;
	}

	//free memory
	free(b);
	free(gam);
	free(Gam);
	free(Phi);
	free(Rho);
	free(Pi);
	free(ps);
	free(nY2);
	free(ps1);
	free(nY21);
	free(Gram2);
	free(ps2);
	gsl_matrix_free(matrix);
	gsl_permutation_free(permutation);
	free(XiPhiR);
	free(YiRhoR);
	free(gam2);
	free(pi2);
	free(X2);
	free(Y2);
	free(sqNorm2);
}
