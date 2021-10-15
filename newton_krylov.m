function [h,converged,k_newton,k_GMRES] = newton_krylov(h_n,F,options)
%NEWTON_KRYLOV Newton-Krylov solver that numerically approximates the
%solution vector h that solves F(h)=0 at the current time-step with a range
%of techniques including GMRES and line searching
% Inputs:
%   h_n: solution vector at previous time-step
%   F: function handle nonlinear system function with solution vector as
%      input
%   options: structure containing structures that define the constants for:
%            - optionsNewton,
%            - optionsLineSearching,
%            - optionsGMRES,
%            - optionsJacobian.
% Outputs:
%   h: approximation of solution vector that solves F(h)=0
%   k: number of newton iterations taken until sufficient convergence

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
newtonStepMethod = optionsNewton.newtonStepMethod;

% Initialise
F_k = F(h_n);
norm_Fk = norm(F_k,2);
norm_F0 = norm_Fk;
h_k = h_n;
k_GMRES = 0;
converged = false;
h = zeros(length(h_n),1);

% Perform Newton iterations until convergence
for k_newton = 1:maxiters
    k_newton
    % Generate Jacobian every m iterations
    if ~mod(k_newton-1,m)
        J = @(h) jacobian(h,F,F_k,optionsJacobian);
    end
    
    % Solve for the Newton step
    if isequal(newtonStepMethod, 'backslash')
        delta_h = -J(h_k)\F_k;
    elseif isequal(newtonStepMethod, 'GMRES')
        [delta_h,k_GMRES,~] = gmres_pre(-J(h_k),F_k,h_k,optionsGMRES);
        if isequal(delta_h,[])
            break;
        end
    else
        error("Unknown newton step generation method");
    end
    
    % Use Line Searching to converge if full newton step insufficient
    if k_newton > 5
        lambda = lineSearching(h_k,delta_h,F,norm_Fk,optionsLineSearching);
    else
        lambda = 1;
    end
    
    if lambda == -1
        break;
    else
        % Accept Newton iteration
        h_k = h_k + lambda*delta_h;
        F_k = F(h_k);
        norm_Fk = norm(F_k,2)
    end
    
    % Check for convergence
    if norm_Fk <= rtol*norm_F0 + atol
        h = h_k;
        converged = true;
        break;
    end
end

end

