% Exercise 4
% --------------------------------------------------------------------------------------

%    1. Create a 2x2 covariance matrix Σ with non-zero off-diagonal elements.
%    2. Generate a white Gaussian random vector w with dimension 2 × 1000, 
%           where each column has a zero mean and an identity covariance matrix.
%    3. Perform the Cholesky decomposition on Σ and obtain its lower triangular matrix L.
%    4. Use the matrix L to map w to a new Gaussian random vector x such that x = Lw.
%    5. Compute the sample covariance matrix of x and verify that the sample covariance 
%       matrix is approximately equal to Σ.


% The covariance matrix is always both symmetric and positive semi- definite.
% Since the off diagonal elements of Σ are fixed to be non zero the
% covariance matrix becomes symmetric and positive definite.


A = [3 7 ; 9 11]
E = (A*A')./2

tf = issymmetric(E)
d = eig(E)
isposdef = all(d > 0)

L = chol(E,'lower')

w = wgn(2,1000,1)
cov_W = cov(w(1,:),w(2,:))

x = L*w

cov_X = cov(x(1,:),x(2,:))
n = norm(E-cov_X,"fro")