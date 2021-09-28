function F = Ffunc(h,h_n,k,psi,Q_p,q_rain,nodes,meshConfig,discretisationConsts)
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
delta = meshConfig.delta;
Delta = meshConfig.Delta;
DV = meshConfig.DV;
quadMats = meshConfig.quadMats;
matNames = meshConfig.matNames;
K = meshConfig.K;

% Extract discretisation constants
dt = discretisationConsts.dt;
theta = discretisationConsts.theta;
Kc = discretisationConsts.Kc;
Hc = discretisationConsts.Hc;
Xc = discretisationConsts.Xc;

% Extract nodes
x = nodes.xNodes;
z = nodes.zNodes;

% Define readability constants
east = 1; west = 2; north = 3; south = 4;
Nx = length(x); Nz = length(z);

% Initialise vectors used in time-discretisation F
G = zeros(Nx*Nz,1); G_n = zeros(Nx*Nz,1);
Psi = zeros(Nx*Nz,1); Psi_n = zeros(Nx*Nz,1);
Q = zeros(Nx*Nz,1); Q_n = zeros(Nx*Nz,1);

% Iterate over the mesh nodes to generate Psi_n, G_n, Q_n, Psi, G, Q
for ii = 1:Nx
    q_east = 0;
    q_west = 0;
    q_north = 0;
    q_south = 0;
    for jj = 1:Nz
        
        % Generate H, psi and k functions for a given h value
        H = @(h) h + z(jj);
        hpMats = squeeze(matNames==quadMats(jj,ii,:)); % materials in surrounding quadrants of node
        psi_h = @(h) sum((psi(h)*hpMats).*(squeeze(DV(jj,ii,:))'))/(Delta(jj,ii,1)*Delta(jj,ii,2));
        k_h = @(h) sum((k(h)*hpMats).*(squeeze(DV(jj,ii,:))'))/(Delta(jj,ii,1)*Delta(jj,ii,2));
        
        % USE h TO GENERATE Psi, G and Q ----------------------------------
        
        % Generate flux values for each face of interior nodes control volumes
        if ii < Nx
            sigma_e = (H(h(jj,ii))<H(h(jj,ii+1))) + (1/2)*(H(h(jj,ii))==H(h(jj,ii+1)));
            k_e = (1-sigma_e)*k_h(h(jj,ii)) + sigma_e*k_h(h(jj,ii+1));
            q_east = k_e * K(jj,ii,east) * (H(h(jj,ii+1)) - H(h(jj,ii)))/delta(jj,ii,east);
        end
        if ii > 1
            sigma_w = (H(h(jj,ii))>H(h(jj,ii-1))) + (1/2)*(H(h(jj,ii))==H(h(jj,ii-1)));
            k_w = (1-sigma_w)*k_h(h(jj,ii-1)) + sigma_w*k_h(h(jj,ii));
            q_west = - k_w * K(jj,ii,west) * (H(h(jj,ii)) - H(h(jj,ii-1)))/delta(jj,ii,west);
        end
        if jj < Nz
            sigma_n = (H(h(jj,ii))<H(h(jj+1,ii))) + (1/2)*(H(h(jj,ii))==H(h(jj+1,ii)));
            k_n = (1-sigma_n)*k_h(h(jj,ii)) + sigma_n*k_h(h(jj+1,ii));
            q_north = k_n * K(jj,ii,north) * (H(h(jj+1,ii)) - H(h(jj,ii)))/delta(jj,ii,north);
        end
        if jj > 1
            sigma_s = (H(h(jj,ii))>H(h(jj-1,ii))) + (1/2)*(H(h(jj,ii))==H(h(jj-1,ii)));
            k_s = (1-sigma_s)*k_h(h(jj-1,ii)) + sigma_s*k_h(h(jj,ii));
            q_south = - k_s * K(jj,ii,south) * (H(h(jj,ii)) - H(h(jj-1,ii)))/delta(jj,ii,south);
        end
        
        % Apply boundary conditions
        q_east = (ii < Nx)*q_east;
        q_west = (ii > 1)*q_west + ...
            (ii == 1 & 3<=z(jj) & z(jj)<=5 & H(h(jj,ii)) > Hc) * ...
            (Kc*(Hc - H(h(jj,ii)))/Xc);
        q_north = (jj < Nz)*q_north + (jj == Nz)*-q_rain;
        q_south = (jj > 1)*q_south;
        
        % Generate and save Q value at current node
        % psi_p(h>=0) generates average psi_sat value in node domain
        Q(Nz*(ii-1)+jj) = (psi_h(h(jj,ii))/psi_h(0) > 0.5) * Q_p(x(ii),z(jj));
        % Generate and save G value at current node
        G(Nz*(ii-1)+jj) = (q_east + q_west + q_north + q_south)/(Delta(jj,ii,1)*Delta(jj,ii,2));
        % Generate and save Psi value at current node
        Psi(Nz*(ii-1)+jj) = psi_h(h(jj,ii));
        
        % USE h_n TO GENERATE Psi_n, G_n and Q_n --------------------------
        
        % Generate flux values for each face of node control volume
        if ii < Nx
            k_e = (1-sigma_e)*k_h(h_n(jj,ii)) + sigma_e*k_h(h_n(jj,ii+1));
            q_east = k_e * K(jj,ii,east) * (H(h_n(jj,ii+1)) - H(h_n(jj,ii)))/delta(jj,ii,east);
        end
        if ii > 1
            k_w = (1-sigma_w)*k_h(h_n(jj,ii-1)) + sigma_w*k_h(h_n(jj,ii));
            q_west = - k_w * K(jj,ii,west) * (H(h_n(jj,ii)) - H(h_n(jj,ii-1)))/delta(jj,ii,west);
        end
        if jj < Nz
            k_n = (1-sigma_n)*k_h(h_n(jj,ii)) + sigma_n*k_h(h_n(jj+1,ii));
            q_north = k_n * K(jj,ii,north) * (H(h_n(jj+1,ii)) - H(h_n(jj,ii)))/delta(jj,ii,north);
        end
        if jj > 1
            k_s = (1-sigma_s)*k_h(h_n(jj-1,ii)) + sigma_s*k_h(h_n(jj,ii));
            q_south = - k_s * K(jj,ii,south) * (H(h_n(jj,ii)) - H(h_n(jj-1,ii)))/delta(jj,ii,south);
        end
        
        % Apply boundary conditions
        q_east = (ii < Nx)*q_east;
        q_west = (ii > 1)*q_west + ...
            (ii == 1 & 3<=z(jj) & z(jj)<=5 & H(h_n(jj,ii)) > Hc) * ...
            (Kc*(Hc - H(h_n(jj,ii)))/Xc);
        q_north = (jj < Nz)*q_north + (jj == Nz)*-q_rain;
        q_south = (jj > 1)*q_south;
        
        % Generate and save Q value at current node
        % psi_p(h>=0) generates average psi_sat value in node domain
        Q_n(Nz*(ii-1)+jj) = (psi_h(h_n(jj,ii))/psi_h(0) > 0.5) * Q_p(x(ii),z(jj));
        % Generate and save G value at current node
        G_n(Nz*(ii-1)+jj) = (q_east + q_west + q_north + q_south)/(Delta(jj,ii,1)*Delta(jj,ii,2));
        % Generate and save Psi value at current node
        Psi_n(Nz*(ii-1)+jj) = psi_h(h_n(jj,ii));
        
    end
end

% Form F with Psi, G, Q, Psi_n, G_n, Q_n
F = Psi - Psi_n - dt*((theta*G + (1-theta)*G_n) + (theta*Q + (1-theta)*Q_n));

end
