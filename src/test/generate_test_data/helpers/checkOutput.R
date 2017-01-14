checkOutput = function(varName, array, refArray, tol)
{
	print(paste("Checking ",varName,sep=""))
	maxError = max(abs(array - refArray))
	if(maxError >= tol)
	{
		print(paste("Inaccuracy: max(abs(error)) = ",maxError," >= ",tol,sep=""))
	} else
	{
		print("OK")
	}
}
