function invXY = CholeskyInverse(X,Y)
%___________________
% CholeskyInverse.m
%
% compute inv(X)*Y for symmetric psd matrix X using cholesky decomposition
% and back substitution
% if X is an array of matrices, do X(:,:,k)\Y for all k in the array
% if Y is empty, return inv(X)

% $Id: CholeskyInverse.m v0.02 2012-06-23 14:45:16 fj $ readjusted to original non-fix version
    
%%%% propaganda
% $$$ myLogo						= cafe_logo( mfilename) ;

num_mats = size(X,3);

if nargin == 1
    Y = eye(size(X(:,:,1)));
end

for k = 1 : num_mats
    L = chol(X(:,:,k), 'lower');
    invXY(:,:,k) = zeros(size(X,2),size(Y,2));
    for j = 1 : size(Y,2)
        v = L\Y(:,j);
        invXY(:,j,k) = L'\v;
    end
end

%%%% Firdaus fix... 
function invXY = fixCholeskyInverse(X,Y)
% temporary bug fix - 2012-06-16 FJ
num_mats = size(X,3);
invXY = zeros(size(X));

if nargin == 1
    Y = eye(size(X(:,:,1)));
end

for k = 1 : num_mats
    if rcond(X(:,:,k)) < eps
        invXY(:,:,k) = X(:,:,k)\Y;
    else
        invXY(:,:,k) = (X(:,:,k)+eps*eye(size(X(:,:,1))))\Y;
    end
end


function invXY = realCholeskyInverse(X,Y)
% compute inv(X)*Y for symmetric psd matrix X using cholesky decomposition
% and back substitution
% if X is an array of matrices, do X(:,:,k)\Y for all k in the array
% if Y is empty, return inv(X)

num_mats = size(X,3);

if nargin == 1
    Y = eye(size(X(:,:,1)));
end

for k = 1 : num_mats
    L = chol(X(:,:,k), 'lower');
    invXY(:,:,k) = zeros(size(X,2),size(Y,2));
    for j = 1 : size(Y,2)
        v = L\Y(:,j);
        invXY(:,j,k) = L'\v;
    end
end


function X=backsub(A,B)

%Input    - A is an n x n upper-triangular nonsingular matrix
%	         - B is an n x 1 matrix
%Output - X is the solution to the linear system AX = B

%  NUMERICAL METHODS: Matlab Programs
% (c) 2004 by John H. Mathews and Kurtis D. Fink
%  Complementary Software to accompany the textbook:
%  NUMERICAL METHODS: Using Matlab, Fourth Edition
%  ISBN: 0-13-065248-2
%  Prentice-Hall Pub. Inc.
%  One Lake Street
%  Upper Saddle River, NJ 07458

%Find the dimension of B and initialize X
 n=length(B);
 X=zeros(n,1);
 X(n)=B(n)/A(n,n);

for k=n-1:-1:1
 X(k)=(B(k)-A(k,k+1:n)*X(k+1:n))/A(k,k);
end