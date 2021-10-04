function [Psi,flux] = FfuncVars_bucket(h,nodes,k,psi,K_vals,quadMats,DV,deltas,Deltas,discretisationConsts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x = nodes.xNodes;
z = nodes.zNodes;
Nx = length(x); Nz = length(z);

q_rain = discretisationConsts.q_rain;

K_e = K_vals.east;
K_w = K_vals.west;
K_n = K_vals.north;
K_s = K_vals.south;

Delta_x = Deltas.x;
Delta_z = Deltas.z;
Delta_xz = Deltas.xz;

delta_e = deltas.east;
delta_w = deltas.west;
delta_n = deltas.north;
delta_s = deltas.south;


% Vectors containing neighbouring heads
h_e = NaN(Nx*Nz,1); h_e(1:end-Nz) = h(Nz+1:end);
h_w = NaN(Nx*Nz,1); h_w(Nz+1:end) = h(1:end-Nz);
h_n = circshift(h,-1); h_n(Nz:Nz:end) = NaN;
h_s = circshift(h,1); h_s(1:Nz:end) = NaN;

% Total head with vectors of neighbouring total heads
H = h + repmat(z,Nx,1);
H_e = NaN(Nx*Nz,1); H_e(1:end-Nz) = H(Nz+1:end);
H_w = NaN(Nx*Nz,1); H_w(Nz+1:end) = H(1:end-Nz);
H_n = circshift(H,-1); H_n(Nz:Nz:end) = NaN;
H_s = circshift(H,1); H_s(1:Nz:end) = NaN;

% Sigma values for each node in each direction
sigma_e = (H<H_e);
sigma_w = (H>H_w);
sigma_n = (H<H_n);
sigma_s = (H>H_s);

% Approximate k at each node for each surrounding material
k_hp = k(h)'; k_he = k(h_e)'; k_hw = k(h_w)'; k_hn = k(h_n)'; k_hs = k(h_s)';
k_hp = sum((k_hp(quadMats).*DV),2) ./ Delta_xz;
k_he = sum((k_he(quadMats).*DV),2) ./ Delta_xz;
k_hw = sum((k_hw(quadMats).*DV),2) ./ Delta_xz;
k_hn = sum((k_hn(quadMats).*DV),2) ./ Delta_xz;
k_hs = sum((k_hs(quadMats).*DV),2) ./ Delta_xz;
k_e = (1-sigma_e).*k_hp + sigma_e.*k_he; k_e(isnan(k_e)) = 0;
k_w = (1-sigma_w).*k_hw + sigma_w.*k_hp; k_w(isnan(k_w)) = 0;
k_n = (1-sigma_n).*k_hp + sigma_n.*k_hn; k_n(isnan(k_n)) = 0;
k_s = (1-sigma_s).*k_hs + sigma_s.*k_hp; k_s(isnan(k_s)) = 0;

% Flux values for each node in each direction
q_e = -k_e .* K_e .* (H_e - H) ./ delta_e;
q_w = k_w .* K_w .* (H - H_w) ./ delta_w;
q_n = -k_n .* K_n .* (H_n - H) ./ delta_n;
q_s = k_s .* K_s .* (H - H_s) ./ delta_s;

% Apply boundary conditions
q_e(isnan(q_e)) = 0;
q_w(isnan(q_w)) = 0;
q_n(isnan(q_n)) = -q_rain;
q_s(isnan(q_s)) = 0;

% Form G, Q and Psi
% need to form a Nx*Nz by 4 (because 4 mats) array containing index of
% material needed at each node for k(h) and psi(h)
flux = (q_e + q_w + q_n + q_s) ./ (Delta_xz);

psi_h = psi(h)';
Psi = sum((psi_h(quadMats).*DV),2) ./ Delta_xz;

end

