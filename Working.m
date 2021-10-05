clear; close all; clc;
% profile on;
% Define material rectangle limits/bounds (x1,x2;z1,z2)
alluvium = [0,5;0,1]; alluvium(:,:,2) = [0,10;1,2]; alluvium(:,:,3) = [0,15;2,10];
materials.alluvium = alluvium;
materials.cappingLayer = [15,100;6.5,8.5];
materials.topsoil = [15,100;8.5,10];
landfill = [5,10;0,1]; landfill(:,:,2) = [10,15;0,2]; landfill(:,:,3) = [15,100;0,6.5];
materials.landfill = landfill;
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
xNodes = sort(unique(xNodes))'; zNodes = sort(unique(zNodes))';
Nx = length(xNodes); Nz = length(zNodes);
matNames = fieldnames(materials); Nmats = length(matNames);

% Plot nodes on material distribution plot
nodesTotal = Nx*Nz;
for x = xNodes'
    for z = zNodes'
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
[K_vals,S,k,psi,Q,deltas,Deltas,DV,quadMats] = nodeConstants2D(materials,constants,xNodes,zNodes);

% Collate in structure for easy transfer between functions
meshConfig.deltas = deltas;
meshConfig.Deltas = Deltas;
meshConfig.DV = DV;
meshConfig.K_vals = K_vals;

% Remove zeros from quadMats so no indexing errors
% These edits will not affect computation as DV is 0 where quadMats is 0
quadMats(quadMats==0) = 1;
idxCorrection = repmat(0:Nmats:Nmats*Nx*Nz-1,Nmats,1)';
quadMats = quadMats + idxCorrection;
meshConfig.quadMats = quadMats;

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
discretisationConsts.q_rain = 0.05; % just a random constant choice

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


%% Solve over time

figure;

h_solved = zeros(Nz,Nx,length(t));
S_solved = zeros(Nz,Nx,length(t));
psi_solved = zeros(Nz,Nx,length(t));

h_solved(:,:,1) = -1 + ((-5 + 1)*repmat(zNodes',Nx,1)')/L2;
h_n = reshape(h_solved(:,:,1),[Nz*Nx,1]);

avgSatsMeasured = zeros(length(t),1);

for t_n = 2:length(t)
    
    % Save psi and S values of previous time-step
    psi_h_n = psi(h_n)';
    psi_h_n = sum((psi_h_n(quadMats).*DV),2) ./ Deltas.xz;
    psi_solved(:,:,t_n-1) = reshape(psi_h_n,[Nz,Nx]);
    S_h_n = S(h_n)';
    S_h_n = sum((S_h_n(quadMats).*DV),2) ./ Deltas.xz;
    S_solved(:,:,t_n-1) = reshape(S_h_n,[Nz,Nx]);
    
    % Plot the previous time solution and store the values as vector
    % heads
    solutionPlot(1) = subplot(2,3,1);
    surf(xNodes,zNodes,h_solved(:,:,t_n-1));
    view(2)
    colormap(solutionPlot(1),flipud(autumn))
    shading interp;
    colorbar
    % water content
    solutionPlot(2) = subplot(2,3,2);
    surf(xNodes,zNodes,psi_solved(:,:,t_n-1));
    view(2)
    colormap(solutionPlot(2),winter)
    shading interp;
    colorbar
    % saturation
    solutionPlot(3) = subplot(2,3,3);
    surf(xNodes,zNodes,S_solved(:,:,t_n-1));
    view(2)
    colormap(solutionPlot(3),cool)
    shading interp;
    colorbar
    % average water content (moisture)
    avgSatsMeasured(t_n-1) = sum(psi_h_n .* Deltas.xz)/(L1*L2);
    solutionPlot(4) = subplot(2,3,[4,5,6]);
    plot(t(1:t_n-1),avgSatsMeasured(1:t_n-1),'b');
    ylim([0,0.5])
    title("average water content at time-step " + num2str(t_n-1));
    drawnow;
    
    % Form F for current time-step
    F = @(h) Ffunc(h,h_n,k,psi,Q,nodes,meshConfig,discretisationConsts);
    
    % Obtain h for current time-step with Newton-Krylov
    [h_n,~] = newton_krylov(h_n,F,options);
    
    % Reshape vector solution and store in solution matrix
    h_solved(:,:,t_n) = reshape(h_n,[Nz,Nx]);
    
end


%%
% close all;
% figure
% 
% for i = 1:40
%     plot(i,i,'r*')
%     hold on;
% plot(1:i,1:i,'k')
% 
% 
% xlim([0 50])
% ylim([0 50])
% drawnow;
% pause(0.25)
% end
% 
% 
% %%
% 
% A = [1,2,3,4;
%     5,6,7,8];
% [Nj,Ni] = size(A);
% A_vec = zeros(Ni*Nj,1);
% for i = 1:Ni
%     for j = 1:Nj
%         A_vec(Nj*(i-1)+j) = A(j,i);
%     end
% end
% A_vec
% 
% 
% 
% 
