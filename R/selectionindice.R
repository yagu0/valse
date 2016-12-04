selectionindice = function(phi, seuil)
{
	dim_phi = dim(phi)
	pp = dim_phi[1]
	m = dim_phi[2]

	A = matrix(0, pp, m)
	B = matrix(0, pp, m)

	for(j in 1:pp)
	{
		cpt1 = 0
		cpt2 = 0
		for(mm in 1:m)
		{
			if(max(phi[j,mm,]) > seuil)
			{
				cpt1 = cpt1 + 1
				A[j,cpt] = mm
			} else
			{
				cpt2 = cpt2+1
				B[j, cpt2] = mm
			}
		}
	}
	return (list(A,B))
}
