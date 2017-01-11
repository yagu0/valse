%covariance matrix for tests on synthetic data: A(i,j) = a^|i-j|
function[A] = covariance(p,a)

	A = a*ones(p,p);
	for i=1:p
		A(i,:) = A(i,:) .^ abs(i-(1:p));
	end

end
