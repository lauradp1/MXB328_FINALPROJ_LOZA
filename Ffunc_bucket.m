function F = Ffunc_bucket(h,h_n,k,psi,nodes,meshConfig,discretisationConsts)
%FFUNC Generates a discretisation of the 2-dimensional PDE of Richard's
%equation at the current time-step to be solved via Newton-Krylov
% Inputs:
%   h: heads at current time-step
%   h_n: heads at previous time-step
%   k: function handle vector of length no. of materials with h as input
%      and outputs the value of k for each material in same order as the
%      keys/names of the materials
%   psi: function handle vector of length no. of materials with h as input
%        and outputs the value of psi for each material in same order as
%        the keys/names of the materials
%   Q: function handle with x and z as input and outputs evapotranspiration
%      term at the given x and z position
%   q_rain: value of the rainfall flux at current time-step
%   nodes: structure containing xNodes and zNodes for the mesh
%   meshConfig: structure containing arrays relating to different aspects
%               of the mesh:
%           - delta: Nx by Nz by 4 array containing values for delta_east,
%                    delta_west, ... at every node
%           - Delta: Nx by Nz by 2 array containing the values for Delta_x
%                    and Delta_z at every node
%           - DV: Nx by Nz by 4 array containing values for the areas of
%                 each quadrant V_SCV surrounding every node
%           - quadMats: Nx by Nz by 4 string array containing the materials
%                       of each quadrant surrounding every node
%           - matNames: cell array containing material names as strings
%           - K: Nx by Nz by 4 array containing the Kxx and Kzz 
%                approximations for each node
%   discretisationConsts: structure containing constants relating to the
%                         discretisation of Richard's PDE:
%           - dt: value of one time-step
%           - theta: 
%           - Kc: hydraulic conductivity of the creek
%           - Hc: total head of the creek
%           - Xc: x-width of the creek
% Outputs:
%   F: Evaluation of F(u) = 0 for the 2-dimensional PDE of Richard's eqn

% Extract meshConfiguration variables
deltas = meshConfig.deltas;
Deltas = meshConfig.Deltas;
DV = meshConfig.DV;
K_vals = meshConfig.K_vals;
quadMats = meshConfig.quadMats;

% Extract discretisation constants
dt = discretisationConsts.dt;
theta = discretisationConsts.theta;

% Obtain Psi, flux and Q for h
[Psi,flux] = FfuncVars_bucket(h,nodes,k,psi,K_vals,quadMats,DV,deltas,Deltas,discretisationConsts);

% Obtain Psi, flux and Q for h_n
[Psi_n,flux_n] = FfuncVars_bucket(h_n,nodes,k,psi,K_vals,quadMats,DV,deltas,Deltas,discretisationConsts);

% Form F with Psi, G, Q, Psi_n, G_n, Q_n
F = Psi - Psi_n + dt*((theta*flux + (1-theta)*flux_n));

end
