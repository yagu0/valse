checkOutput = function(varName, matrix, refMatrix, tol){
  print('Checking %s\n',varName);
  maxError = max(max(max(max(abs(matrix - refMatrix)))));
  if(maxError >= tol){
    print('Inaccuracy: max(abs(error)) = %g >= %g\n',maxError,tol)
  }
  else{
    print('OK\n')
  }
}