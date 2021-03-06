function [h,k] = newton_krylov(h_n,F,options)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Extract options for each function
optionsNewton = options.Newton;
optionsLineSearching = options.LineSearching;
optionsGMRES = options.GMRES;
optionsJacobian = options.Jacobian;

% Extract options for Newton iteration
m = optionsNewton.m;
atol = optionsNewton.atol;
rtol = optionsNewton.rtol;
maxiters = optionsNewton.maxiters;

% Initialise
N = length(h_n);
F_k = F(h_n);
norm_Fk = norm(F_k,2);
norm_F0 = norm_Fk;
h = h_n;

% Perform Newton iterations until convergence
for k = 1:maxiters
    % Generate Jacobian every m iterations
    if ~mod(k-1,m)
        % if jacobian free is just an approximation of J then make it so it
        % gets approximated in the same function here:
        J = jacobian(h,F_k,N,optionsJacobian);
    end
    
    % Solve for the Newton step
        % will J always be in the right format? Or will diagonal need its
        % own variation?
        % is h the right choice for x0?
    delta_h = GMRES(J(h),-F(h),h,optionsGMRES);
    
    % Use Line Searching to converge
    lambda = lineSearching(h,delta_h,F,norm_Fk,optionsLineSearching);
    
    % Accept Newton iteration
    F_k = F(h + lambda*delta_h);
    norm_Fk = norm(F_k,2);
    
    % Check for convergence
    if norm_Fk <= rtol*norm_F0 + atol
        break;
    end
end

end

