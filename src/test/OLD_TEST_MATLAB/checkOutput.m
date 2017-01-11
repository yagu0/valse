function[]=checkOutput(varName, matrix, refMatrix, tol)

	fprintf('Checking %s\n',varName);
	maxError = max(max(max(max(abs(matrix - refMatrix)))));
	if maxError >= tol
		fprintf('    Inaccuracy: max(abs(error)) = %g >= %g\n',maxError,tol);
	else
		fprintf('    OK\n');
	end

end
