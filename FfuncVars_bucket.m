function [Psi,flux] = FfuncVars_bucket(h,nodes,k,psi,K_vals,quadMats,DVs,deltas,Deltas,discretisationConsts)
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

Delta_x_n = Deltas.x_n;
Delta_xz = Deltas.xz;
Delta_xz_e = Deltas.xz_e;
Delta_xz_w = Deltas.xz_w;
Delta_xz_n = Deltas.xz_n;
Delta_xz_s = Deltas.xz_s;

delta_e = deltas.east;
delta_w = deltas.west;
delta_n = deltas.north;
delta_s = deltas.south;

% Vectors containing neighbouring heads
h_e = circshift(h,-Nz); h_e(end-Nz+1:end) = NaN;
h_w = circshift(h,Nz); h_w(1:Nz) = NaN;
h_n = circshift(h,-1); h_n(Nz:Nz:end) = NaN;
h_s = circshift(h,1); h_s(1:Nz:end) = NaN;

% Total head with vectors of neighbouring total heads
H = h + repmat(z,Nx,1);
H_e = circshift(H,-Nz); H_e(end-Nz+1:end) = NaN;
H_w = circshift(H,Nz); H_w(1:Nz) = NaN;
H_n = circshift(H,-1); H_n(Nz:Nz:end) = NaN;
H_s = circshift(H,1); H_s(1:Nz:end) = NaN;

% Generate shifted arrays for DV
DV = DVs.DV;
DV_e = DVs.DV_e;
DV_w = DVs.DV_w;
DV_n = DVs.DV_n;
DV_s = DVs.DV_s;

quadMats_p = quadMats.quadMats_p;
quadMats_e = quadMats.quadMats_e;
quadMats_w = quadMats.quadMats_w;
quadMats_n = quadMats.quadMats_n;
quadMats_s = quadMats.quadMats_s;

% Sigma values for each node in each direction
sigma_e = (H<H_e);
sigma_w = (H>H_w);
sigma_n = (H<H_n);
sigma_s = (H>H_s);

% Approximate k at each node for each surrounding material
k_hp = k(h)'; k_he = k(h_e)'; k_hw = k(h_w)'; k_hn = k(h_n)'; k_hs = k(h_s)';
k_hp = sum((k_hp(quadMats_p).*DV),2) ./ Delta_xz;
k_he = sum((k_he(quadMats_e).*DV_e),2) ./ Delta_xz_e;
k_hw = sum((k_hw(quadMats_w).*DV_w),2) ./ Delta_xz_w;
k_hn = sum((k_hn(quadMats_n).*DV_n),2) ./ Delta_xz_n;
k_hs = sum((k_hs(quadMats_s).*DV_s),2) ./ Delta_xz_s;
k_e = (1-sigma_e).*k_hp + sigma_e.*k_he;
k_w = (1-sigma_w).*k_hw + sigma_w.*k_hp;
k_n = (1-sigma_n).*k_hp + sigma_n.*k_hn;
k_s = (1-sigma_s).*k_hs + sigma_s.*k_hp;

% Flux values for each node in each direction
q_e = -k_e .* K_e .* (H_e - H) ./ delta_e;
q_w = k_w .* K_w .* (H - H_w) ./ delta_w;
q_n = -k_n .* K_n .* (H_n - H) ./ delta_n;
q_s = k_s .* K_s .* (H - H_s) ./ delta_s;

% Apply boundary conditions
q_e(isnan(q_e)) = 0;
q_w(isnan(q_w)) = 0;
q_n(isnan(q_n)) = -q_rain * Delta_x_n(isnan(q_n));
q_s(isnan(q_s)) = 0;

% Form G, Q and Psi
% need to form a Nx*Nz by 4 (because 4 mats) array containing index of
% material needed at each node for k(h) and psi(h)
flux = (q_e + q_w + q_n + q_s) ./ (Delta_xz);

psi_h = psi(h)';
Psi = sum((psi_h(quadMats_p).*DV),2) ./ Delta_xz;

end

