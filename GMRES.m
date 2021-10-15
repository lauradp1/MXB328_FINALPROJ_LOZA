function [x,m] = GMRES(A,b,x0,options)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Extract options
atol = options.atol;
rtol = options.rtol;
maxiters = options.maxiters;
precond = options.precond;

% Handle precondition
if isequal(precond,'Jacobi')
    M = diag(diag(A));
elseif isequal(precond,'Gauss-Seidel')
    M = tril(A);
end

% Initialise
N = size(A,1);
H = zeros(maxiters+1, maxiters);
V = zeros(N, maxiters+1);
rnorm = zeros(maxiters, 1);

% Right-preconditioned FOM
r = b - A*x0; beta = norm(r,2); V(:,1) = r/beta; converged = false; rnorm0 = beta;
e = @(m) [beta; zeros(m,1)];

for m = 1:maxiters
    V(:,m+1) = A*(M\V(:,m));
    for j = 1:m
        H(j,m) = V(:,j)'*V(:,m+1);
        V(:,m+1) = V(:,m+1) - H(j,m)*V(:,j);
    end
    H(m+1,m) = norm(V(:,m+1),2);
    
    if H(m+1,m) < 1e-14
        y = H(1:m,1:m) \ e(m-1);
        converged = true;
        break;
    else
        V(:,m+1) = V(:,m+1)/H(m+1,m);
    end
    
    y = H(1:(m+1),1:m) \ e(m);
    rnorm(m) = norm(e(m) - H(1:m+1,1:m)*y,2);
    
    if rnorm(m) < rtol*rnorm0 + atol
        converged = true;
        break;
    end
end

if converged
    x = x0 + M\(V(:,1:m)*y);
else
    x = [];
end

end

