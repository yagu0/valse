#include "EMGLLF.h"
#include <gsl/gsl_linalg.h>

// TODO: don't recompute indexes every time......
void EMGLLF(
	// IN parameters
	const double* phiInit,  // parametre initial de moyenne renormalisé
	const double* rhoInit,  // parametre initial de variance renormalisé
	const double* piInit,   // parametre initial des proportions
	const double* gamInit,  // paramètre initial des probabilités a posteriori de chaque échantillon
	int mini,      // nombre minimal d'itérations dans l'algorithme EM
	int maxi,      // nombre maximal d'itérations dans l'algorithme EM
	double gamma,  // valeur de gamma : puissance des proportions dans la pénalisation pour un Lasso adaptatif
	double lambda, // valeur du paramètre de régularisation du Lasso
	const double* X,     // régresseurs
	const double* Y,     // réponse
	double tau,    // seuil pour accepter la convergence
	// OUT parameters (all pointers, to be modified)
	double* phi,  // parametre de moyenne renormalisé, calculé par l'EM
	double* rho,  // parametre de variance renormalisé, calculé par l'EM
	double* pi,   // parametre des proportions renormalisé, calculé par l'EM
	double* LLF,   // log vraisemblance associé à cet échantillon, pour les valeurs estimées des paramètres
    double* S,
	// additional size parameters
	int n,         // nombre d'echantillons
	int p,         // nombre de covariables
	int m,         // taille de Y (multivarié)
	int k)         // nombre de composantes dans le mélange
{
	//Initialize outputs
	copyArray(phiInit, phi, p*m*k);
	copyArray(rhoInit, rho, m*m*k);
	copyArray(piInit, pi, k);
	zeroArray(LLF, maxi);
	//S is already allocated, and doesn't need to be 'zeroed'

	//Other local variables
	//NOTE: variables order is always [maxi],n,p,m,k
	double* gam = (double*)malloc(n*k*sizeof(double));
	copyArray(gamInit, gam, n*k);
	double* b = (double*)malloc(k*sizeof(double));
	double* Phi = (double*)malloc(p*m*k*sizeof(double));
	double* Rho = (double*)malloc(m*m*k*sizeof(double));
	double* Pi = (double*)malloc(k*sizeof(double));
	double* gam2 = (double*)malloc(k*sizeof(double));
	double* pi2 = (double*)malloc(k*sizeof(double));
	double* Gram2 = (double*)malloc(p*p*k*sizeof(double));
	double* ps = (double*)malloc(m*k*sizeof(double));
	double* nY2 = (double*)malloc(m*k*sizeof(double));
	double* ps1 = (double*)malloc(n*m*k*sizeof(double));
	double* ps2 = (double*)malloc(p*m*k*sizeof(double));
	double* nY21 = (double*)malloc(n*m*k*sizeof(double));
	double* Gam = (double*)malloc(n*k*sizeof(double));
	double* X2 = (double*)malloc(n*p*k*sizeof(double));
	double* Y2 = (double*)malloc(n*m*k*sizeof(double));
	gsl_matrix* matrix = gsl_matrix_alloc(m, m);
	gsl_permutation* permutation = gsl_permutation_alloc(m);
	double* YiRhoR = (double*)malloc(m*sizeof(double));
	double* XiPhiR = (double*)malloc(m*sizeof(double));
	double dist = 0.;
	double dist2 = 0.;
	int ite = 0;
	double EPS = 1e-15;
	double* dotProducts = (double*)malloc(k*sizeof(double));

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
				//Y2(:,mm,r)=sqrt(gam(:,r)).*transpose(Y(mm,:));
				for (int u=0; u<n; u++)
					Y2[ai(u,mm,r,n,m,k)] = sqrt(gam[mi(u,r,n,k)]) * Y[mi(u,mm,m,n)];
			}
			for (int i=0; i<n; i++)
			{
				//X2(i,:,r)=X(i,:).*sqrt(gam(i,r));
				for (int u=0; u<p; u++)
					X2[ai(i,u,r,n,m,k)] = sqrt(gam[mi(i,r,n,k)]) * X[mi(i,u,n,p)];
			}
			for (int mm=0; mm<m; mm++)
			{
				//ps2(:,mm,r)=transpose(X2(:,:,r))*Y2(:,mm,r);
				for (int u=0; u<p; u++)
				{
					double dotProduct = 0.;
					for (int v=0; v<n; v++)
						dotProduct += X2[ai(v,u,r,n,m,k)] * Y2[ai(v,mm,r,n,m,k)];
					ps2[ai(u,mm,r,n,m,k)] = dotProduct;
				}
			}
			for (int j=0; j<p; j++)
			{
				for (int s=0; s<p; s++)
				{
					//Gram2(j,s,r)=transpose(X2(:,j,r))*(X2(:,s,r));
					double dotProduct = 0.;
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
			//b(r) = sum(sum(abs(phi(:,:,r))));
			double sumAbsPhi = 0.;
			for (int u=0; u<p; u++)
				for (int v=0; v<m; v++)
					sumAbsPhi += fabs(phi[ai(u,v,r,p,m,k)]);
			b[r] = sumAbsPhi;
		}
		//gam2 = sum(gam,1);
		for (int u=0; u<k; u++)
		{
			double sumOnColumn = 0.;
			for (int v=0; v<n; v++)
				sumOnColumn += gam[mi(v,u,n,k)];
			gam2[u] = sumOnColumn;
		}
		//a=sum(gam*transpose(log(pi)));
		double a = 0.;
		for (int u=0; u<n; u++)
		{
			double dotProduct = 0.;
			for (int v=0; v<k; v++)
				dotProduct += gam[mi(u,v,n,k)] * log(pi[v]);
			a += dotProduct;
		}

		//tant que les proportions sont negatives
		int kk = 0;
		int pi2AllPositive = 0;
		double invN = 1./n;
		while (!pi2AllPositive)
		{
			//pi2(:)=pi(:)+0.1^kk*(1/n*gam2(:)-pi(:));
			for (int r=0; r<k; r++)
				pi2[r] = pi[r] + pow(0.1,kk) * (invN*gam2[r] - pi[r]);
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

		//t(m) la plus grande valeur dans la grille O.1^k tel que ce soit décroissante ou constante
		//(pi.^gamma)*b
		double piPowGammaDotB = 0.;
		for (int v=0; v<k; v++)
			piPowGammaDotB += pow(pi[v],gamma) * b[v];
		//(pi2.^gamma)*b
		double pi2PowGammaDotB = 0.;
		for (int v=0; v<k; v++)
			pi2PowGammaDotB += pow(pi2[v],gamma) * b[v];
		//transpose(gam2)*log(pi2)
		double prodGam2logPi2 = 0.;
		for (int v=0; v<k; v++)
			prodGam2logPi2 += gam2[v] * log(pi2[v]);
		while (-invN*a + lambda*piPowGammaDotB < -invN*prodGam2logPi2 + lambda*pi2PowGammaDotB && kk<1000)
		{
			//pi2=pi+0.1^kk*(1/n*gam2-pi);
			for (int v=0; v<k; v++)
				pi2[v] = pi[v] + pow(0.1,kk) * (invN*gam2[v] - pi[v]);
			//pi2 was updated, so we recompute pi2PowGammaDotB and prodGam2logPi2
			pi2PowGammaDotB = 0.;
			for (int v=0; v<k; v++)
				pi2PowGammaDotB += pow(pi2[v],gamma) * b[v];
			prodGam2logPi2 = 0.;
			for (int v=0; v<k; v++)
				prodGam2logPi2 += gam2[v] * log(pi2[v]);
			kk++;
		}
		double t = pow(0.1,kk);
		//sum(pi+t*(pi2-pi))
		double sumPiPlusTbyDiff = 0.;
		for (int v=0; v<k; v++)
			sumPiPlusTbyDiff += (pi[v] + t*(pi2[v] - pi[v]));
		//pi=(pi+t*(pi2-pi))/sum(pi+t*(pi2-pi));
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
					double dotProduct = 0.0;
					for (int u=0; u<p; u++)
						dotProduct += X2[ai(i,u,r,n,p,k)] * phi[ai(u,mm,r,n,m,k)];
					//ps1(i,mm,r)=Y2(i,mm,r)*dot(X2(i,:,r),phi(:,mm,r));
					ps1[ai(i,mm,r,n,m,k)] = Y2[ai(i,mm,r,n,m,k)] * dotProduct;
					nY21[ai(i,mm,r,n,m,k)] = Y2[ai(i,mm,r,n,m,k)] * Y2[ai(i,mm,r,n,m,k)];
				}
				//ps(mm,r)=sum(ps1(:,mm,r));
				double sumPs1 = 0.0;
				for (int u=0; u<n; u++)
					sumPs1 += ps1[ai(u,mm,r,n,m,k)];
				ps[mi(mm,r,m,k)] = sumPs1;
				//nY2(mm,r)=sum(nY21(:,mm,r));
				double sumNy21 = 0.0;
				for (int u=0; u<n; u++)
					sumNy21 += nY21[ai(u,mm,r,n,m,k)];
				nY2[mi(mm,r,m,k)] = sumNy21;
				//rho(mm,mm,r)=((ps(mm,r)+sqrt(ps(mm,r)^2+4*nY2(mm,r)*(gam2(r))))/(2*nY2(mm,r)));
				rho[ai(mm,mm,k,m,m,k)] = ( ps[mi(mm,r,m,k)] + sqrt( ps[mi(mm,r,m,k)]*ps[mi(mm,r,m,k)] 
					+ 4*nY2[mi(mm,r,m,k)] * (gam2[r]) ) ) / (2*nY2[mi(mm,r,m,k)]);
			}
		}
		for (int r=0; r<k; r++)
		{
			for (int j=0; j<p; j++)
			{
				for (int mm=0; mm<m; mm++)
				{
					//sum(phi(1:j-1,mm,r).*transpose(Gram2(j,1:j-1,r)))+sum(phi(j+1:p,mm,r).*transpose(Gram2(j,j+1:p,r)))
					double dotPhiGram2 = 0.0;
					for (int u=0; u<j; u++)
						dotPhiGram2 += phi[ai(u,mm,r,p,m,k)] * Gram2[ai(j,u,r,p,p,k)];
					for (int u=j+1; u<p; u++)
						dotPhiGram2 += phi[ai(u,mm,r,p,m,k)] * Gram2[ai(j,u,r,p,p,k)];
					//S(j,r,mm)=-rho(mm,mm,r)*ps2(j,mm,r)+sum(phi(1:j-1,mm,r).*transpose(Gram2(j,1:j-1,r)))
					//    +sum(phi(j+1:p,mm,r).*transpose(Gram2(j,j+1:p,r)));
					S[ai(j,mm,r,p,m,k)] = -rho[ai(mm,mm,r,m,m,k)] * ps2[ai(j,mm,r,p,m,k)] + dotPhiGram2;
					if (fabs(S[ai(j,mm,r,p,m,k)]) <= n*lambda*pow(pi[r],gamma))
						phi[ai(j,mm,r,p,m,k)] = 0;
					else if (S[ai(j,mm,r,p,m,k)] > n*lambda*pow(pi[r],gamma))
						phi[ai(j,mm,r,p,m,k)] = (n*lambda*pow(pi[r],gamma) - S[ai(j,mm,r,p,m,k)]) 
							/ Gram2[ai(j,j,r,p,p,k)];
					else
						phi[ai(j,mm,r,p,m,k)] = -(n*lambda*pow(pi[r],gamma) + S[ai(j,mm,r,p,m,k)]) 
							/ Gram2[ai(j,j,r,p,p,k)];
				}
			}
		}

		/////////////
		// Etape E //
		/////////////

		int signum;
		double sumLogLLF2 = 0.0;
		for (int i=0; i<n; i++)
		{
			double sumLLF1 = 0.0;
			double sumGamI = 0.0;
			double minDotProduct = INFINITY;

			for (int r=0; r<k; r++)
			{
				//Compute
				//Gam(i,r) = Pi(r) * det(Rho(:,:,r)) * exp( -1/2 * (Y(i,:)*Rho(:,:,r) - X(i,:)...
				//    *phi(:,:,r)) * transpose( Y(i,:)*Rho(:,:,r) - X(i,:)*phi(:,:,r) ) );
				//split in several sub-steps
				
				//compute Y(i,:)*rho(:,:,r)
				for (int u=0; u<m; u++)
				{
					YiRhoR[u] = 0.0;
					for (int v=0; v<m; v++)
						YiRhoR[u] += Y[imi(i,v,n,m)] * rho[ai(v,u,r,m,m,k)];
				}

				//compute X(i,:)*phi(:,:,r)
				for (int u=0; u<m; u++)
				{
					XiPhiR[u] = 0.0;
					for (int v=0; v<p; v++)
						XiPhiR[u] += X[mi(i,v,n,p)] * phi[ai(v,u,r,p,m,k)];
				}

				// compute dotProduct < Y(:,i)*rho(:,:,r)-X(i,:)*phi(:,:,r) . Y(:,i)*rho(:,:,r)-X(i,:)*phi(:,:,r) >
				dotProducts[r] = 0.0;
				for (int u=0; u<m; u++)
					dotProducts[r] += (YiRhoR[u]-XiPhiR[u]) * (YiRhoR[u]-XiPhiR[u]);
				if (dotProducts[r] < minDotProduct)
					minDotProduct = dotProducts[r];
			}
			double shift = 0.5*minDotProduct;
			for (int r=0; r<k; r++)
			{
				//compute det(rho(:,:,r)) [TODO: avoid re-computations]
				for (int u=0; u<m; u++)
				{
					for (int v=0; v<m; v++)
						matrix->data[u*m+v] = rho[ai(u,v,r,m,m,k)];
				}
				gsl_linalg_LU_decomp(matrix, permutation, &signum);
				double detRhoR = gsl_linalg_LU_det(matrix, signum);

				Gam[mi(i,r,n,k)] = pi[r] * detRhoR * exp(-0.5*dotProducts[r] + shift);
				sumLLF1 += Gam[mi(i,r,n,k)] / pow(2*M_PI,m/2.0);
				sumGamI += Gam[mi(i,r,n,k)];
			}
			sumLogLLF2 += log(sumLLF1);
			for (int r=0; r<k; r++)
			{
				//gam(i,r)=Gam(i,r)/sum(Gam(i,:));
				gam[mi(i,r,n,k)] = sumGamI > EPS
					? Gam[mi(i,r,n,k)] / sumGamI
					: 0.0;
			}
		}
		
		//sum(pen(ite,:))
		double sumPen = 0.0;
		for (int r=0; r<k; r++)
			sumPen += pow(pi[r],gamma) * b[r];
		//LLF(ite)=-1/n*sum(log(LLF2(ite,:)))+lambda*sum(pen(ite,:));
		LLF[ite] = -invN * sumLogLLF2 + lambda * sumPen;
		if (ite == 0)
			dist = LLF[ite];
		else 
			dist = (LLF[ite] - LLF[ite-1]) / (1.0 + fabs(LLF[ite]));
		
		//Dist1=max(max((abs(phi-Phi))./(1+abs(phi))));
		double Dist1 = 0.0;
		for (int u=0; u<p; u++)
		{
			for (int v=0; v<m; v++)
			{
				for (int w=0; w<k; w++)
				{
					double tmpDist = fabs(phi[ai(u,v,w,p,m,k)]-Phi[ai(u,v,w,p,m,k)]) 
						/ (1.0+fabs(phi[ai(u,v,w,p,m,k)]));
					if (tmpDist > Dist1)
						Dist1 = tmpDist;
				}
			}
		}
		//Dist2=max(max((abs(rho-Rho))./(1+abs(rho))));
		double Dist2 = 0.0;
		for (int u=0; u<m; u++)
		{
			for (int v=0; v<m; v++)
			{
				for (int w=0; w<k; w++)
				{
					double tmpDist = fabs(rho[ai(u,v,w,m,m,k)]-Rho[ai(u,v,w,m,m,k)]) 
						/ (1.0+fabs(rho[ai(u,v,w,m,m,k)]));
					if (tmpDist > Dist2)
						Dist2 = tmpDist;
				}
			}
		}
		//Dist3=max(max((abs(pi-Pi))./(1+abs(Pi))));
		double Dist3 = 0.0;
		for (int u=0; u<n; u++)
		{
			for (int v=0; v<k; v++)
			{
				double tmpDist = fabs(pi[v]-Pi[v]) / (1.0+fabs(pi[v]));
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
	free(dotProducts);
}
