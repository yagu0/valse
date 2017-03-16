#' Among a collection of models, this function constructs a subcollection of models with
#' models having strictly different dimensions, keeping the model which minimizes
#' the likelihood if there were several with the same dimension
#'
#' @param LLF a matrix, the first column corresponds to likelihoods for several models
#'				the second column corresponds to the dimensions of the corresponding models.
#'
#' @return a list with indices, a vector of indices selected models,
#'				 and D1, a vector of corresponding dimensions
#' @export
#'
modelSelection = function(LLF)
{
	D = LLF[,2]
	D1 = unique(D)

	indices = rep(1, length(D1))
	#select argmax MLE
	if (length(D1)>2)
	{
		for (i in 1:length(D1))
		{
			A = c()
			for (j in 1:length(D))
			{
				if(D[[j]]==D1[[i]])
					a = c(a, LLF[j,1])
			}
			b = max(a)
			#indices[i] : first indices of the binary vector where u_i ==1
			indices[i] = which.max(LLF == b)
		}
	}

	return (list(indices=indices,D1=D1))
}

#TODO:
## Programme qui sélectionne un modèle
## proposer à l'utilisation différents critères (BIC, AIC, slope heuristic)
