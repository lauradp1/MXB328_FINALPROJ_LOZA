function [K,S,k,psi,Q,delta,Delta,DV,quadMats] = nodeConstants2D(materials,constants,xNodes,zNodes)
%NODECONSTANTS2D Computes time-independent variables for given mesh so to
%provide fast access while iterating over time
% Inputs:
%   materials: structure containing 2 by 2 by n arrays with each matrix
%              consisting of the x limits (1,:) and z limits (2,:) of the
%              2D material for n 2D boxes of the material
%   constants: structure with keys for each constant where each key maps to
%              a vector of length the number of materials containing the
%              different values of the constant for each material
%   xNodes: vector of x-values corresponding to the x-nodes
%   zNodes: vector of z-values corresponding to the z-nodes
% Outputs:
%   K: Nx by Nz by 4 array containing the Kxx and Kzz approximations for
%      each node
%   k: function handle vector of length no. of materials with h as input
%      and outputs the value of k for each material in same order as the
%      keys/names of the materials
%   psi: function handle vector of length no. of materials with h as input
%        and outputs the value of psi for each material in same order as
%        the keys/names of the materials
%   Q: function handle with x and z as input and outputs evapotranspiration
%      term at the given x and z position
%   delta: Nx by Nz by 2 by 2 array containing distances from each node to
%          the surrounding nodes
%   Delta: Nx by Nz by 2 array containing the length and height of each
%          node domain Vp
%   DV: Nx by Nz by 4 array containing areas of the 4 quadrants surrounding
%       each node
%   quadMats: Nx by Nz by 4 array containing the materials of the quadrants
%             surrounding each node

% Define readability constants
east = 1; west = 2; north = 3; south = 4;
Nx = length(xNodes); Nz = length(zNodes);

% Extract material constants
Kxx = constants.Kxx;
Kzz = constants.Kzz;
psi_res = constants.psi_res;
psi_sat = constants.psi_sat;
alpha = constants.alpha;
n = constants.n;
m = constants.m;
l = constants.l; l1 = l(1); l2 = l(2);
R = constants.R; R1 = R(1); R2 = R(2);
L = constants.L; L1 = L(1); L2 = L(2);
matNames = fieldnames(materials); % Material names

% Generate S(h), k(h), psi(h) and Q(x,z) for each material
S = @(h) (h>=0) + (h<0)*(1 + (-alpha*h).^n).^(-m);
k = @(h) (h>=0) + (h<0)*sqrt(S(h)).*(1 - (1 - S(h).^(m.^(-1))).^m ).^2;
psi = @(h) (h>=0)*psi_sat + (h<0)*(psi_res + S(h).*(psi_sat - psi_res));
Q = @(x,z) (0<=x & x<=15 & (L2-l1)<=z & z<=L2)*((L2-l1<=z & z<=L2)*(-R1*(z-L2+l1)^2)/l1^2) + ...
    (15<x & x<=L1 & (L2-l2)<=z & z<=L2)*((L2-l2<=z & z<=L2)*(-R2*(z-L2+l2)^2)/l2^2);

% Iterate through each node and calculate time-independant values e.g. Kxx
% and Kzz approximations to use at each node East, West, North and South
K = zeros(Nz,Nx,4); % 4:Kxx_east,Kxx_west,Kzz_north,Kzz_south
delta = zeros(Nz,Nx,4); % 4:dx_east,dx_west,dz_north,dz_south
Delta = zeros(Nz,Nx,2); % 2:Deltax_p,Deltaz_p
DV = zeros(Nz,Nx,4); % 4:DV_quadrant1,DV_quadrant2,...
quadMats = strings(Nz,Nx,4);
for i = 1:Nx
    x_i = xNodes(i);    % Save x-value of current node
    % Save East and West node values
    if i < Nx; x_e = xNodes(i+1); else; x_e = x_i; end
    if i > 1; x_w = xNodes(i-1); else; x_w = x_i; end
    for j = 1:Nz
        z_i = zNodes(j);    % Save z-value of current node
        % Save North and South node values
        if j > 1; z_s = zNodes(j-1); else; z_s = z_i; end
        if j < Nz; z_n = zNodes(j+1); else; z_n = z_i; end
        
        % Compute deltax and deltaz values for East, West, North and South
        delta(j,i,:) = [x_e-x_i,x_i-x_w,z_n-z_i,z_i-z_s];
        % Compute Deltax and Deltaz
        Delta(j,i,:) = [(delta(j,i,east)+delta(j,i,west))/2,(delta(j,i,north)+delta(j,i,south))/2];
        
        % Determine Kxx and Kzz for each quadrant (find material of each
        % corner of the current node domain Vp
        nodeMats = zeros(2,4);
        if x_e ~= x_i && z_n ~= z_i     % quadrant 1
            nodeMat = convertCharsToStrings(nodeMaterial(materials,...
                [x_i+delta(j,i,east)/2,z_i+delta(j,i,north)/2]));
            nodeMats(:,1) = [Kxx(matNames==nodeMat) Kzz(matNames==nodeMat)];
            % Variables for k_p(h) and psi_p(h) approximations
            DV(j,i,1) = (delta(j,i,east)*delta(j,i,north))/4;
            quadMats(j,i,1) = nodeMat;
        end
        if x_w ~= x_i && z_n ~= z_i     % quadrant 2
            nodeMat = convertCharsToStrings(nodeMaterial(materials,...
                [x_i-delta(j,i,west)/2,z_i+delta(j,i,north)/2]));
            nodeMats(:,2) = [Kxx(matNames==nodeMat) Kzz(matNames==nodeMat)];
            % Variables for k_p(h) and psi_p(h) approximations
            DV(j,i,2) = (delta(j,i,west)*delta(j,i,north))/4;
            quadMats(j,i,2) = nodeMat;
        end
        if x_w ~= x_i && z_s ~= z_i     % quadrant 3
            nodeMat = convertCharsToStrings(nodeMaterial(materials,...
                [x_i-delta(j,i,west)/2,z_i-delta(j,i,south)/2]));
            nodeMats(:,3) = [Kxx(matNames==nodeMat) Kzz(matNames==nodeMat)];
            % Variables for k_p(h) and psi_p(h) approximations
            DV(j,i,3) = (delta(j,i,west)*delta(j,i,south))/4;
            quadMats(j,i,3) = nodeMat;
        end
        if x_e ~= x_i && z_s ~= z_i     % quadrant 4
            nodeMat = convertCharsToStrings(nodeMaterial(materials,...
                [x_i+delta(j,i,east)/2,z_i-delta(j,i,south)/2]));
            nodeMats(:,4) = [Kxx(matNames==nodeMat) Kzz(matNames==nodeMat)];
            % Variables for k_p(h) and psi_p(h) approximations
            DV(j,i,4) = (delta(j,i,east)*delta(j,i,south))/4;
            quadMats(j,i,4) = nodeMat;
        end
        
        % Compute the Kxx (East,West) and Kzz (North,South) approximations
        % THIS MAY BE WRONG SHOULD DOUBLE CHECK VALUES ARE WHAT WE WANT
        K(j,i,east) = (-nodeMats(1,4)*delta(j,i,south) - ...
            nodeMats(1,1)*delta(j,i,north)) / (2);
        K(j,i,west) = (-nodeMats(1,3)*delta(j,i,south) - ...
            nodeMats(1,2)*delta(j,i,north)) / (2);
        K(j,i,north) = (-nodeMats(2,1)*delta(j,i,east) - ...
            nodeMats(2,2)*delta(j,i,west)) / (2);
        K(j,i,south) = (-nodeMats(2,4)*delta(j,i,east) - ...
            nodeMats(2,3)*delta(j,i,west)) / (2);
    end
end

end

