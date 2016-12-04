suppressionmodelesegaux2 = function(B1,rho,pi)
{
	ind = c()
	dim_B1 = dim(B1)
	B2 = array(0,dim=c(dim_B1[1],dim_B1[2],dim_B1[3]))
	nombreLambda=dim_B1[[2]]
	glambda = rep(0,nombreLambda)

	#for(j in 1:nombreLambda){
	#	for(ll in 1:(l-1)){
	#		if(B1[,,l] == B1[,,ll]){
	#			ind = c(ind, l)
	#		}
	#	}
	#}
	#ind = unique(ind)
	#B1 = B1[,,-ind]
	#rho = rho[,,,-ind] 
	#pi = pi[,-ind]

	suppressmodel = suppressionmodelesegaux(B1,B2,glambda,rho,pi)
	return (list(B1 = suppressmodel$B1, ind = suppressmodel$B2,
		rho = suppressmodel$rho, pi = suppressmodel$pi))
}
