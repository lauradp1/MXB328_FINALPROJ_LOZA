function [F,outflow] = Ffunc(h,h_n,k,psi,Q_p,nodes,meshConfig,discretisationConsts)
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
%   Q_p: function handle with q_rain constant as input and outputs vector
%        corresponding to evapotranspiration term at each node
%   nodes: structure containing xNodes and zNodes for the mesh
%   meshConfig: structure containing arrays relating to different aspects
%               of the mesh:
%           - deltas: structure containing delta vectors for east, west,
%                     north and south nodes
%           - Deltas: structure containing Delta vectors for x and z nodes
%                     along with Delta_x.*Delta_z for nodes and east, west,
%                     north and south nodes. Also Delta_creek for the creek
%                     boundary condition
%           - DVs: structure containing DV and DV shifted for east, west,
%                  north and south nodes
%           - K_vals: structure containing K  for east, west, north and
%                     south nodes
%           - quadMats: structure containing indexing arrays for which
%                       materials are at which nodes along with the shifted
%                       versions for east, west, north and south nodes
%   discretisationConsts: structure containing constants relating to the
%                         discretisation of Richard's PDE:
%           - dt: value of time-step
%           - theta: theta method constant
%           - Kc: hydraulic conductivity of the creek
%           - Hc: total head of the creek
%           - Xc: x-width of the creek
%           - q_rain: value of the rainfall flux at current time-step
% Outputs:
%   F: Evaluation of F(u) = 0 for the 2-dimensional PDE of Richard's eqn

% Extract meshConfiguration variables
deltas = meshConfig.deltas;
Deltas = meshConfig.Deltas;
DVs = meshConfig.DVs;
K_vals = meshConfig.K_vals;
quadMats = meshConfig.quadMats;

% Extract discretisation constants
dt = discretisationConsts.dt;
theta = discretisationConsts.theta;

% Obtain Psi, flux and Q for h
[Psi,flux,Q,outflow] = FfuncVars(h,nodes,k,psi,Q_p,K_vals,quadMats,DVs,deltas,Deltas,discretisationConsts);

% Obtain Psi, flux and Q for h_n
[Psi_n,flux_n,Q_n,~] = FfuncVars(h_n,nodes,k,psi,Q_p,K_vals,quadMats,DVs,deltas,Deltas,discretisationConsts);

% Form F with Psi, G, Q, Psi_n, G_n, Q_n
F = Psi - Psi_n + dt*((theta*flux + (1-theta)*flux_n) - (theta*Q + (1-theta)*Q_n));

end
