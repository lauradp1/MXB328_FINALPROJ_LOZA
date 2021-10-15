function [Psi,flux,Q] = FfuncVars(h,nodes,k,psi,Q_p,K_vals,quadMats,DVs,deltas,Deltas,discretisationConsts)
%FFUNCVARS Computes Psi, flux and Q vectors for use with Ffunc
% Inputs:
%   h: Nx*Nz column vector of solution heads
%   nodes: structure containing xNodes and zNodes vectors
%   k: function handle vector of length no. of materials with h as input
%      and outputs the value of k for each material in same order as the
%      keys/names of the materials
%   psi: function handle vector of length no. of materials with h as input
%        and outputs the value of psi for each material in same order as
%        the keys/names of the materials
%   Q_p: function handle with q_rain constant as input and outputs vector
%        corresponding to evapotranspiration term at each node
%   K_vals: structure containing K vectors for east, west, north and south
%           nodes corresponding to the hydraulic conductivity at each node
%   quadMats: structure containing indexing arrays for nodes and east, 
%             west, north and south shifted nodes which define which
%             materials surround each node
%   DVs: structure containing DV (areas of surrounding control volumes) for
%        nodes and east, west, north and south shifted nodes
%   deltas: structure containing delta vectors for east, west, north and
%           south shifted nodes (used for k approximations)
%   Deltas: structure containing Delta vectors for x and z dimension along
%           with Delta_x*Delta_z and it's east, west, north and south
%           shifts and Delta_creek for the west boundary condition
%   discretisationConsts: structure containing boundary condition
%                         constants:
%           - Kc: hydraulic conductivity of the creek
%           - Hc: total head of the creek
%           - Xc: x-width of the creek
%           - q_rain: value of the rainfall flux at current time-step
% Outputs:
%   Psi: Nx*Nz column vector containing the psi value at each node
%   flux: Nx*Nz column vector containing the flux value at each node
%   Q: Nx*Nz column vector containing the Q value at each node

% Extract nodes
x = nodes.xNodes;
z = nodes.zNodes;
Nx = length(x); Nz = length(z);

% Extract boundary condition constants
Kc = discretisationConsts.Kc;
Hc = discretisationConsts.Hc;
Xc = discretisationConsts.Xc;
q_rain = discretisationConsts.q_rain;

% Extract K vectors
K_e = K_vals.east;
K_w = K_vals.west;
K_n = K_vals.north;
K_s = K_vals.south;

% Extract Delta vectors
Delta_x = Deltas.x;
Delta_creek = Deltas.creek;
Delta_xz = Deltas.xz;
Delta_xz_e = Deltas.xz_e;
Delta_xz_w = Deltas.xz_w;
Delta_xz_n = Deltas.xz_n;
Delta_xz_s = Deltas.xz_s;

% Extract delta vectors
delta_e = deltas.east;
delta_w = deltas.west;
delta_n = deltas.north;
delta_s = deltas.south;

% Extract DV arrays
DV = DVs.DV;
DV_e = DVs.DV_e;
DV_w = DVs.DV_w;
DV_n = DVs.DV_n;
DV_s = DVs.DV_s;

% Extract quadMats arrays
quadMats_p = quadMats.quadMats_p;
quadMats_e = quadMats.quadMats_e;
quadMats_w = quadMats.quadMats_w;
quadMats_n = quadMats.quadMats_n;
quadMats_s = quadMats.quadMats_s;

% head vectors containing neighbouring heads
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
q_w(isnan(q_w)) = (3<=z & z<=5 & H(1:Nz)>Hc) .* (Kc*(Hc - H(1:Nz))/Xc) ...
                    .* Delta_creek(1:Nz);
q_n(isnan(q_n)) = -q_rain * Delta_x(isnan(q_n));
q_s(isnan(q_s)) = 0;

% Form flux, Psi and Q
flux = (q_e + q_w + q_n + q_s) ./ (Delta_xz);

psi_h = psi(h)';
Psi = sum((psi_h(quadMats_p).*DV),2) ./ Delta_xz;

psi_sat = psi(zeros(Nx*Nz,1))';
Psi_sat = sum((psi_sat(quadMats_p).*DV),2) ./ Delta_xz;
Q = ((Psi./Psi_sat) > 0.5) .* Q_p(q_rain);

end

