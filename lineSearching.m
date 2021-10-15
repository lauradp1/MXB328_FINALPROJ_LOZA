function lambda = lineSearching(h,delta_h,F,norm_Fk,options)
%LINESEARCHING 2-point parabolic linesearching computes a lambda such that
%h + lambda*delta_h is sufficiently convergent
% Inputs:
%   h: solution vector for previous newton step
%   delta_h: newton step vector
%   F: function handle nonlinear system function that takes a h vector as
%      input
%   norm_Fk: norm of F(h) computed for this newton iteration
%   options: structure containing linesearching constants:
%            - dev: linesearching tolerance to prevent lambda tending to
%                   zero
% Outputs:
%   lambda: constant used to define the approximate solution vector h = h +
%           lambda*delta_h

% Extract options and initialise norm
dev = options.dev;
lambda = 1;
F_k = F(h + lambda*delta_h);
norm_ls = norm(F_k,2);

while norm_ls >= (1 - dev*lambda)*norm_Fk
    g0 = norm_Fk^2;
    glambda = norm_ls^2;
    lambda_new = g0*lambda^2/(glambda + (2*lambda - 1)*g0); 
    
    if lambda_new < 0.1*lambda
        lambda = 0.1*lambda;
    elseif lambda_new > 0.5*lambda
        lambda = 0.5*lambda;
    else
        lambda = lambda_new;
    end
    
    h_k = h + lambda*delta_h;
    F_k = F(h_k);
    norm_ls = norm(F_k,2);
    
    if lambda < dev
        lambda = -1;
        break;
    end
end

end


