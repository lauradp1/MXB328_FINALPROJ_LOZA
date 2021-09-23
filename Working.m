clear; close all; clc;
% profile on;
% Define material rectangle limits/bounds (x1,x2;z1,z2)
alluvium = [0,5;0,1]; alluvium(:,:,2) = [0,10;1,2]; alluvium(:,:,3) = [0,15;2,10];
materials.alluvium = alluvium;
materials.cappingLayer = [15,100;6.5,8.5];
materials.topsoil = [15,100;8.5,10];
landfill = [5,10;0,1]; landfill(:,:,2) = [10,15;0,2]; landfill(:,:,3) = [15,100;0,6.5];
materials.landfill = landfill;
matNames = fieldnames(materials); % Save names of materials for later
% Note: algorithm does not check for gaps in material distribution

% Define mesh parameters [x,z]
c_r = [1,0.2];        % Refinement constants
r = [1.2,1];          % Constants to determine spread of nodes (after refinement)
numNodes = [8,20];    % Number of x and z nodes to add (after refinement)

% Generate base mesh
[nodes,refinements,meshLims,materialsPlot] = meshMaterialNodes(materials,c_r);
xNodes = nodes{1}; zNodes = nodes{2};
L1 = meshLims(2,1);
L2 = meshLims(2,2);

% Fill mesh with some more nodes (xMin replaced with 20 as 0 to 20 is
% already sufficiently populated with nodes)
meshLims_temp = meshLims; meshLims_temp(1,1) = 20;
[nodes_fill] = meshNodes(meshLims_temp,r,numNodes,refinements);
xNodes = [xNodes nodes_fill{1}]; zNodes = [zNodes nodes_fill{2}];
xNodes = sort(unique(xNodes)); zNodes = sort(unique(zNodes));

% Plot nodes on material distribution plot
nodesTotal = length(xNodes)*length(zNodes);
for x = xNodes
    for z = zNodes
        plot(x, z, 'k.');
        hold on
    end
end
set(materialsPlot, 'visible', 'on');
% profview
%%
% clear; clc;
% Define constants in same order as materials:
%   alluvium, cappingLayer, topsoil, landfill
constants.Kxx = [17.28,0.052,0.2,10];
constants.Kzz = [0.4*17.28,0.4*0.052,0.2*0.2,0.2*10];
constants.psi_res = [0.01,0.2,0.05,0.025];
constants.psi_sat = [0.33,0.47,0.15,0.38];
constants.alpha = [1.43,1.04,0.95,1.5];
constants.n = [1.51,1.3954,1.256,1.2];
constants.m = [1-1/1.51,1-1/1.3954,1-1/1.256,1-1/1.2];
constants.l = [1.5,3];
constants.R = [0.1,0.2];
constants.L = [L1,L2];

% Compute time-independent variables
[K,k,psi,Q,delta,Delta,DV,quadMats] = nodeConstants2D(materials,constants,xNodes,zNodes);

% Collate in structure for easy transfer between functions
meshConfig.delta = delta;
meshConfig.Delta = Delta;
meshConfig.DV = DV;
meshConfig.quadMats = quadMats;
meshConfig.matNames = matNames;
meshConfig.K = K;

% Collate nodes
clear nodes;
nodes.xNodes = xNodes;
nodes.zNodes = zNodes;

% Collate discretisation constants
discretisationConsts.dt = 1; % should replace this with = dt; and define dt at top with t
discretisationConsts.theta = [];
discretisationConsts.Kc = [];
discretisationConsts.Hc = [];
discretisationConsts.Xc = [];

q_rain = 4; % just a random constant choice
h_0 = -1 + ((-5 + 1).*repmat(zNodes,length(xNodes),1))./L2;

F = @(h) Ffunc(h,h_0,k,psi,Q,q_rain,nodes,meshConfig,discretisationConsts);


%%

% Nmats by 4 to tell what material at each quadrant
% hpMats = squeeze(matNames==quadMats(10,10,:));
% 
% h = -2;
% 
% % This is to compute k_hp and psi_hp approximations
% psi_p = @(h) sum((psi(h)*hpMats).*(squeeze(DV(10,30,:))'))/(Delta(10,30,1)*Delta(10,30,2));
% k_p = sum((k(h)*hpMats).*(squeeze(DV(10,30,:))'))/(Delta(10,30,1)*Delta(10,30,2));
% 
% % 10,10 would be replaced with i,j in loop in Gfunc
% 
% % For using Q will need to also do (psi_p(h(i,j))/psi_p(0) > 0.5)*Q
% % since psi_p(0) (or any h >= 0) corresponds to only psi_sat averages)
% 
% 
% (psi_p(h)/psi_p(0) > 0.5)*Q(x(i),z(i))












