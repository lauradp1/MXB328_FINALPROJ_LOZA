function lambda = lineSearching(h,delta_h,F,norm,options)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Extract options and initialise norm
dev = options.dev;
lambda = 1;
F_k = F(h + lambda*delta_h);
norm_Fk = norm(F_k,2);

while norm_Fk >= (1 - dev*lambda)*norm
    g0 = norm^2;
    glambda = norm_Fk^2;
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
    norm_Fk = norm(F_k,2);
    
    if lambda < dev
        error('Line searching aborted as lambda is less than dev');
    end
end

end

