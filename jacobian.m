function J = jacobian(h_n,F,F_k,options)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[Nz,Nx] = size(h_n);
N = Nx*Nz;
J = zeros(N,N);     % Initialise J
% Variables
    norm_hn = norm(h_n,2);
    
    % Approximate J
    for j = 1:N
        % determine epsilon
        if norm_hn == 0
            epsilon = sqrt(eps);
        else
            epsilon = sqrt(eps)*norm_hn;
        end
        % determine J at current j
        h_nCopy = h_n;
        h_nCopy(j) = h_nCopy(j) + epsilon;
        J(:,j) = (F(h_nCopy) - F_k)/epsilon;
    end

end

