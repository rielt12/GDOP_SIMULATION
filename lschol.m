function x = lschol(A,b)
%LSCHOL   LEAST SQUARES USING CHOLESKY DECOMPOSITION.
%   x = LSCHOL(A,b) computes the n-dimensional column vector x that 
%   minimizes norm(b-A*x), where A is an m-by-n coefficient matrix and 
%   b is the m-dimensional right side column vector ( m >> n ).
%
%   See also CHOL, LSQR, MLDIVIDE
%
% % Example:
%
% m = 1e3;
% n = 3;
% A = rand(m,n); % A is a tall matrix
% b = rand(m,1);
% disp('Computing least squares using Cholesky decomposition');
% tic; xchol = lschol(A,b); toc;
% disp('Now computing least squares using QR decomposition (just for comparison)');
% tic; xqr = lsqr(A,b); toc;
% xchol
% xqr

%   Copyright Marco Cococcioni, 2016

% Compute the Cholesky decomposition
L = chol(A'*A); % L is upper triangular and (A'A == L'*L)
d = A'*b;

n = length(d);
% Solve L'z = d (via backwards substitution on a lower triangular matrix)
z = backward_subtitution_assuming_lower_triangular(L',d,n);

% Solve Lx=z (via backwards substitution on an upper triangular matrix)
x = backward_subtitution_assuming_upper_triangular(L,z,n);

end


function x = backward_subtitution_assuming_upper_triangular(T,y,n)
x=zeros(n,1);
x(n)=y(n)/T(n,n);
for i=n-1:-1:1
    x(i)=(y(i)-T(i,i+1:n)*x(i+1:n))/T(i,i);
end
end

function x = backward_subtitution_assuming_lower_triangular(T,y,n)
x=zeros(n,1);
x(1)=y(1)/T(1,1);
for i=2:n
    x(i)=(y(i)-T(i,1:i-1)*x(1:i-1))/T(i,i);
end
end
