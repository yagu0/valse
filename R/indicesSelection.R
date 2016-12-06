#' Construct the set of relevant indices -> ED: je crois que cette fonction n'est pas utile
#'
#' @param phi regression matrix, of size p*m
#' @param thresh threshold to say a cofficient is equal to zero
#'
#' @return a list with A, a matrix with relevant indices (size = p*m) and B, a 
#'          matrix with irrelevant indices (size = p*m)
#' @export
indicesSelection = function(phi, thresh = 1e-6)
{
  dim_phi = dim(phi)
  p = dim_phi[1]
  m = dim_phi[2]
  
  A = matrix(0, p, m)
  B = matrix(0, p, m)
  
  for(j in 1:p)
  {
    cpt1 = 0
    cpt2 = 0
    for(mm in 1:m)
    {
      if(max(phi[j,mm,]) > thresh)
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
  return (list(A=A,B=B))
}
