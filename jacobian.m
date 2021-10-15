function J = jacobian(h_n,F,F_k,options)
%JACOBIAN Computes the full Jacobian for the nonlinear system function F
% Inputs:
%   h_n: Nx*Nz column vector containing solution at previous newton-step
%   F: function handle nonlinear system function that takes a h vector as
%      input
%   F_k: Nx*Nz column vector of already computed F(h_n)
%   options: structure containing constants for Jacobian generation:
%            - Nx: number of x nodes
%            - Nz: number of z nodes
%            - method: string for which Jacobian generation technique to
%                      use (column-wise, banded, jacobian-free)
% Outputs:
%   J: full Jacobian for F of size Nx*Nz by Nx*Nz

% Extract variables and initialise ouput
method = options.method;
Nx = options.Nx;
Nz = options.Nz;
N = Nx*Nz;
J = zeros(N,N);

if isequal(method, 'column-wise')
    % Determine epsilon from h_n
    norm_hn = norm(h_n,2);
    if norm_hn == 0
        epsilon = sqrt(eps);
    else
        epsilon = sqrt(eps)*norm_hn;
    end
    % Compute the Jacobian column-wise
    I = eye(N);
    for j = 1:N
        J(:,j) = (F(h_n + epsilon*I(:,j)) - F_k)/epsilon;
    end
    
elseif isequal(method, 'banded')
    % Compute bandwidth and initialise shift vectors for each band
    bandwidth = 2*Nz + 1;
    shift = zeros(bandwidth,N);
    for band = 1:bandwidth
        shift(band,band:bandwidth:N) = 1;
    end
    % Generate Jacobian bands
    norm_hn = norm(h_n,2);
    J_bands = zeros(bandwidth,N);
    for band = 1:bandwidth
        % Determine epsilon
        if norm_hn == 0
            epsilon = sqrt(eps)/norm(shift(band,:),2);
        else
            epsilon = (sqrt(eps)*norm_hn)/norm(shift(band,:),2);
        end
        % Compute the Jacobian band
        J_bands(band,:) = (F(h_n + epsilon*shift(band,:)') - F_k)/epsilon;
    end
    % Place the bands into the full Jacobian for output
    for i = 1:((bandwidth+1)/2)-1
        J(1:i+((bandwidth+1)/2)-1,i) = J_bands(i,1:i+((bandwidth+1)/2)-1);
    end
    for i = ((bandwidth+1)/2):N - ((bandwidth+1)/2)+1
        if i<=bandwidth
            J((i - (bandwidth+1)/2)+1:(i - (bandwidth+1)/2)+bandwidth, i) = ...
                J_bands(i, (i - (bandwidth+1)/2)+1:(i - (bandwidth+1)/2)+bandwidth)';
        else
            J((i - (bandwidth+1)/2)+1:(i - (bandwidth+1)/2)+bandwidth, i) = ...
                J_bands(mod(i, bandwidth)+(mod(i, bandwidth)==0)*bandwidth, ...
                (i - (bandwidth+1)/2)+1:(i - (bandwidth+1)/2)+bandwidth)';
        end
    end
    for i = N - ((bandwidth+1)/2)+2:N
        J((i - (bandwidth+1)/2)+1:N, i) = J_bands(mod(i, bandwidth) + ...
            (mod(i, bandwidth)==0)*bandwidth, (i - (bandwidth+1)/2)+1:N)';
    end
    
elseif isequal(method, 'jacobian-free')
    
    
    
else
    error("Unknown Jacobian generation method");
end

end

