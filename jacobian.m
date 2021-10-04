function J = jacobian(h_n,F,F_k,options)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[Nz,Nx] = size(h_n);
N = Nx*Nz;
J = zeros(N,N);     % Initialise J
I = eye(N);
% Variables
norm_hn = norm(h_n,2);

if norm_hn == 0
    epsilon = sqrt(eps);
else
    epsilon = sqrt(eps)*norm_hn;
end

% Approximate J
for j = 1:N
    % determine J at current j
    J(:,j) = (F(h_n + epsilon*I(:,j)) - F_k)/epsilon;
end

end

