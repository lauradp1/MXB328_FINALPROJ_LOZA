function [x,m,converged] = gmres_pre(A,b,x0,options)
% Inputs
%   A: coefficent matrix (matrix)
%   b: right hand side vector (vector)
%   x0: initial estimate to solution of Ax = b (vector)
%   atol: absolute error tolerence (positive scalar)
%   rtol: relative error tolerence (positive scalar)
%   maxiters: maximum number of iterations (positive integer)
%   pretype: preconditioner type 'left' or 'right' (string)
%   precond: preconditioner option 'Jacobi' or 'Gauss-Seidel' (string)
% Outputs
%   x: approximate solution of Ax = b (vector) (x = [] is returned if not converged)
%   m: number of iterations taken (positive integer)
%   converged: true if stopping criterion satisfied and false otherwise (logical scalar)

% atol = options.atol;
rtol = options.rtol;
maxiters = options.maxiters;
precond = options.precond;
% pretype = options.pretype;

% Initialise
N  = size(A,1);
H = zeros(maxiters+1,maxiters);
V = zeros(N,maxiters+1);

% GMRES
if isequal(precond, 'Jacobi')
    D = diag(A);
    M = diag(D);
end
if isequal(precond, 'Gauss-Seidel')
    M = tril(A);
    M = full(M);
end
setup.type = 'crout';
setup.milu = 'row';
setup.droptol = 0.03;
[L,U,P]= ilu(sparse(A),setup);

r0= b-A*x0;
% if isequal(pretype, 'left')
%     r0_tilde = M\r0;
%     norm_r0 = norm(r0,2);
%     beta = norm(r0_tilde,2);
%     V(:,1) = r0_tilde/beta;
% else
norm_r0 = norm(r0,2);
beta = norm(r0,2);
V(:,1) = r0/beta; %first basis vector (following gram-schmidt)
resnormp = zeros(maxiters,1);
resnormp(1) = beta;
m =1;
converged = false;
while resnormp(m) > rtol*resnormp(1) && m<maxiters
    v = U\(L\(P*V(:,m)));
    V(:,m+1) = A*v;
% for m = 1:maxiters  
%     if isequal(pretype,'right')
%         V(:,m+1) = A*(M\V(:,m));
%     elseif isequal(pretype, 'left')
%         V(:,m+1) = M\(A*V(:,m));
%     end
    for j = 1:m
        H(j,m) = V(:,j)'*V(:,m+1);
        V(:,m+1)= V(:,m+1) - H(j,m)*V(:,j); %compute next vector for bases V
    end
    H(m+1,m) = norm(V(:,m+1),2); %update hessenberg matrix
    %^this is all the arnoldi sequence
    if H(m+1,m) < 10^-13 %almost zero
        RHS = [beta;zeros(m-1, 1)];
        [U2,S,K] = svd(H(1:m,1:m));
        S = diag(S);
        invS = 1./S;
        invS = diag(invS);
        ym = K*invS*U2'*[beta;zeros(m-1, 1)];
        %ym = H(1:m,1:m)\[beta;zeros(m-1, 1)];
        resnormp(m+1) = norm(RHS-H(1:m,1:m)*ym);
        m = m+1;
        break;
    else
        V(:,m+1) = V(:,m+1)/H(m+1,m);
        RHS = [beta;zeros(m, 1)];
        [U2,S,K] = svd(H(1:m+1,1:m),'econ');
        S = diag(S);
        invS = 1./S;
        invS = diag(invS);
        ym = K*invS*U2'*[beta;zeros(m, 1)];
        %ym = H(1:m+1,1:m)\[beta;zeros(m, 1)];
        resnormp(m+1) = norm(RHS-H(1:m+1,1:m)*ym);
    end
    m = m+1;
end

    
%     %two following lines different for GMRES
%     ym = H(1:m+1, 1:m)\[beta;zeros(m, 1)]; 
%     if isequal(pretype,'left')
%         b= M*beta*V(:,1:m+1)*[1;zeros(m,1)];
%         c = M*V(:,1:m+1)*H(1:m+1, 1:m)*ym;
%         rm =b-c; 
%         norm_rm = norm(rm,2);
%     else
%         rm = [beta;zeros(m,1)] -H(1:m+1, 1:m)*ym; 
%         norm_rm = norm(rm,2);
%     end
%     if norm_rm < rtol*norm_r0 + atol
%         converged = true;
%         break;
%     end
% end

v = U\(L\(P*V(:,1:m-1)*ym));
x = x0+v;
end
% if converged == 0
%     x = [];
% elseif isequal(pretype,'right')
%     x = x0 + M\(V(:, 1:m) * ym);
% else
%     x = x0 + V(:, 1:m) * ym;  
% end