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
c_r = [0,0];
r = [1,1];
numNodes = [5,5];

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
zNodes = zNodes(length(zNodes):-1:1);

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
%% Construct and define constants

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
[K,S,k,psi,Q,delta,Delta,DV,quadMats] = nodeConstants2D(materials,constants,xNodes,zNodes);

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

% Define time range
t = linspace(0,5,60);
dt = t(2)-t(1);

% Collate discretisation constants
discretisationConsts.dt = dt;
discretisationConsts.theta = 0.5;
discretisationConsts.Kc = 0.0108;
discretisationConsts.Hc = 4;
discretisationConsts.Xc = 5;

% Collate Newton method constants
optionsNewton.m = 1;
optionsNewton.atol = 1e-10;
optionsNewton.rtol = 1e-10;
optionsNewton.maxiters = 300;

% Collate Line Searching constants
optionsLineSearching.dev = 1e-4;

% Collate GMRES constants
optionsGMRES.atol = 1e-10;
optionsGMRES.rtol = 1e-10;
optionsGMRES.maxiters = 300;
optionsGMRES.precond = 'Jacobi'; % Jacobi or Gauss-Seidel

% Collate Jacobian constants
optionsJacobian.var = 0;

% Collate all options and constants for Newton-Krylov
options.Newton = optionsNewton;
options.LineSearching = optionsLineSearching;
options.GMRES = optionsGMRES;
options.Jacobian = optionsJacobian;

q_rain = 4; % just a random constant choice


%% Solve over time

figure;
Nx = length(xNodes); Nz = length(zNodes);

h_solved = zeros(Nz,Nx,length(t));
S_solved = zeros(Nz,Nx,length(t));
psi_solved = zeros(Nz,Nx,length(t));

h_solved(:,:,1) = -1 + ((-5 + 1).*repmat(zNodes,length(xNodes),1)')./L2;

for t_n = 2:length(t)
    
    % Plot the previous time solution and store the values as vector
    h_vec = zeros(Nx*Nz,1);
    for i = 1:Nx
        for j = 1:Nz
            hpMats = squeeze(matNames==quadMats(j,i,:)); % materials in surrounding quadrants of node
            S_solved(j,i,t_n-1) = sum((S(h_solved(j,i,t_n-1))*hpMats).*(squeeze(DV(j,i,:))'))/(Delta(j,i,1)*Delta(j,i,2));
            psi_solved(j,i,t_n-1) = sum((psi(h_solved(j,i,t_n-1))*hpMats).*(squeeze(DV(j,i,:))'))/(Delta(j,i,1)*Delta(j,i,2));
            h_vec(Nz*(i-1)+j) = h_solved(j,i,t_n-1);
        end
    end
    % heads
    solutionPlot(1) = subplot(1,3,1);
    surf(xNodes,zNodes,h_solved(:,:,t_n-1));
    view(2)
    colormap(solutionPlot(1),flipud(autumn))
    shading interp;
    colorbar
    % water content
    solutionPlot(2) = subplot(1,3,2);
    surf(xNodes,zNodes,psi_solved(:,:,t_n-1));
    view(2)
    colormap(solutionPlot(2),winter)
    shading interp;
    colorbar
    % saturation
    solutionPlot(3) = subplot(1,3,3);
    surf(xNodes,zNodes,S_solved(:,:,t_n-1));
    view(2)
    colormap(solutionPlot(3),cool)
    shading interp;
    colorbar
    % average water content (moisture)
    avgSat = sum(sum(psi_solved(:,:,t_n-1)))/(Nx*Nz);
    drawnow;
    
    % Form F for current time-step
    h_vec = h_solved(:,:,t_n-1);
    F = @(h) Ffunc(h,h_vec,k,psi,Q,q_rain,nodes,meshConfig,discretisationConsts);
    
    % Obtain h for current time-step with Newton-Krylov
    [h_solved(:,:,t_n),~] = newton_krylov(h_vec,F,options);
    
end


%%
close all;
figure

for i = 1:40
    plot(i,i,'r*')
    hold on;
plot(1:i,1:i,'k')


xlim([0 50])
ylim([0 50])
drawnow;
pause(0.25)
end


%%

A = [1,2,3,4;
    5,6,7,8];
[Nj,Ni] = size(A);
A_vec = zeros(Ni*Nj,1);
for i = 1:Ni
    for j = 1:Nj
        A_vec(Nj*(i-1)+j) = A(j,i);
    end
end
A_vec




