%clearvars -except cosineRain; close all; clc;
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
N = Nx*Nz;
% for x = xNodes'
%     for z = zNodes'
%         plot(x, z, 'k.');
%         hold on
%     end
% end
% set(materialsPlot, 'visible', 'on');
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
constants.creekLims = [3,5];

% Compute time-independent variables
[K_vals,S,k,psi,Q,deltas,Deltas,DVs,quadMats_p] = nodeConstants2D(materials,constants,xNodes,zNodes);

% Collate in structure for easy transfer between functions
meshConfig.deltas = deltas;
meshConfig.Deltas = Deltas;
meshConfig.DVs = DVs;
meshConfig.K_vals = K_vals;

% Remove zeros from quadMats so no indexing errors
% These edits will not affect computation as DV is 0 where quadMats is 0
quadMats_p(quadMats_p==0) = 1;
idxCorrection = repmat(0:Nmats:Nmats*N-1,Nmats,1)';
quadMats.quadMats_p = quadMats_p + idxCorrection;
quadMats.quadMats_e = circshift(quadMats_p,-Nz) + idxCorrection;
quadMats.quadMats_w = circshift(quadMats_p,Nz) + idxCorrection;
quadMats.quadMats_n = circshift(quadMats_p,-1) + idxCorrection;
quadMats.quadMats_s = circshift(quadMats_p,1) + idxCorrection;
quadMats_p = quadMats_p + idxCorrection;
meshConfig.quadMats = quadMats;

% Collate nodes
clear nodes;
nodes.xNodes = xNodes;
nodes.zNodes = zNodes;

% Define time range
t = linspace(0,100,400);
dt = t(2)-t(1);

% Collate discretisation constants
discretisationConsts.dt = dt;
discretisationConsts.theta = 1;
discretisationConsts.Kc = 0.0108;
discretisationConsts.Hc = 4;
discretisationConsts.Xc = 5;
discretisationConsts.L1 = L1;
discretisationConsts.L2 = L2;
discretisationConsts.q_rain = 0.0042;
psi_sat = psi(zeros(Nx*Nz,1))';
Psi_sat = sum((psi_sat(quadMats_p).*DVs.DV),2) ./ Deltas.xz;
discretisationConsts.Psi_sat = Psi_sat;
Psi_avg_max = sum(Psi_sat .* Deltas.xz)/(L1*L2);
discretisationConsts.Psi_avg_max = Psi_avg_max;

% Collate Newton method constants
optionsNewton.m = 1;
optionsNewton.atol = 1e-6;
optionsNewton.rtol = 1e-6;
optionsNewton.maxiters = 30;
optionsNewton.newtonStepMethod = 'backslash'; % GMRES, backslash

% Collate Line Searching constants
optionsLineSearching.dev = 1e-6;

% Collate GMRES constants
optionsGMRES.atol = 1e-8;
optionsGMRES.rtol = 1e-8;
optionsGMRES.maxiters = 50;
optionsGMRES.precond = 'Jacobi'; % Jacobi or Gauss-Seidel
optionsGMRES.pretype = 'right';

% Collate Jacobian constants
optionsJacobian.Nx = Nx;
optionsJacobian.Nz = Nz;
optionsJacobian.method = 'banded'; % banded, column-wise, jacobian-free

% Collate all options and constants for Newton-Krylov
options.Newton = optionsNewton;
options.LineSearching = optionsLineSearching;
options.GMRES = optionsGMRES;
options.Jacobian = optionsJacobian;


%% Solve over time

figure;

t_max = 365*5;
dt = 1;
dt_growthFactor = 1.5;
dt_max = 4;
dt_min = 0.1;
dt_consec = 0;
t = 0; % initialise t
t_vals = 0;
Ksc = 5;
Msc = 30;
Tsc = 5;

h_solved = -1 + ((-5 + 1)*repmat(zNodes',Nx,1)')/L2;
h_n = reshape(h_solved,[N,1]);
psi_h_n = psi(h_n)';
psi_h_n = sum((psi_h_n(quadMats_p).*DVs.DV),2) ./ Deltas.xz;
psi_solved = reshape(psi_h_n,[Nz,Nx]);
S_h_n = S(h_n)';
S_h_n = sum((S_h_n(quadMats_p).*DVs.DV),2) ./ Deltas.xz;
S_solved = reshape(S_h_n,[Nz,Nx]);
avgSatsMeasured = sum(psi_h_n .* Deltas.xz)/(L1*L2);
outflows = zeros(sum(3<=zNodes & zNodes<=5),1);

%%

% Plot solutions at current time-step
solutionPlot(1) = subplot(2,3,1);
contourf(xNodes,zNodes,h_solved(:,:,end));
view(2)
colormap(solutionPlot(1),flipud(autumn))
shading interp;
colorbar
% water content
solutionPlot(2) = subplot(2,3,2);
contourf(xNodes,zNodes,psi_solved(:,:,end));
view(2)
colormap(solutionPlot(2),flipud(winter))
shading interp;
colorbar
% saturation
solutionPlot(3) = subplot(2,3,3);
contourf(xNodes,zNodes,S_solved(:,:,end));
view(2)
colormap(solutionPlot(3),cool)
shading interp;
colorbar
% average water content (moisture)
solutionPlot(4) = subplot(2,3,[4,5,6]);
plot(t_vals,avgSatsMeasured,'b');
ylim([0,0.5])
title("average water content at time " + num2str(t_vals(end)) + " (dt = " + num2str(dt) + ")");
drawnow;

while t + dt < t_max
    converged = false;
    % keep trying until a dt results in a successful convergence
    while ~converged
        discretisationConsts.dt = dt;
        discretisationConsts.q_rain = cosineRain(t + dt);
        discretisationConsts.rainfall = avgSatsMeasured(end)/Psi_avg_max < 0.95;
        discretisationConsts.evapotranspiration = psi_h_n./Psi_sat > 0.5;
        % Form F for current time-step
        F = @(h) Ffunc(h,h_n,k,psi,Q,nodes,meshConfig,discretisationConsts);
        % Obtain h for current time-step with Newton-Krylov
        [h_ntemp,newtonConverged,k_newton,k_GMRES] = newton_krylov(h_n,F,options);
        if isequal(optionsNewton.newtonStepMethod,'backslash'); k_GMRES = Msc; end
        % Determine if dt should be changed
        if ~newtonConverged
            if dt/2 <= dt_min
                error("Minimum dt reached");
            else
                dt = dt/2;
                dt_consec = 0;
            end
        elseif k_newton <= Ksc && k_GMRES <= Msc
            dt_consec = dt_consec + 1;
            if dt*dt_growthFactor < dt_max && dt_consec == Tsc
                dt = dt*dt_growthFactor;
                dt_consec = 0;
            end
            converged = true;
            h_n = h_ntemp;
        else
            converged = true;
            h_n = h_ntemp;
        end
    end
    t_vals(end+1) = t + dt;
    t = t_vals(end)
    
    % Save h, psi and S values of current time-step
    h_solved(:,:,end+1) = reshape(h_n,[Nz,Nx]);
    psi_h_n = psi(h_n)';
    psi_h_n = sum((psi_h_n(quadMats_p).*DVs.DV),2) ./ Deltas.xz;
    psi_solved(:,:,end+1) = reshape(psi_h_n,[Nz,Nx]);
    S_h_n = S(h_n)';
    S_h_n = sum((S_h_n(quadMats_p).*DVs.DV),2) ./ Deltas.xz;
    S_solved(:,:,end+1) = reshape(S_h_n,[Nz,Nx]);
    % save outflow values
    [~,outflows(:,length(t_vals)-1)] = Ffunc(h_n,h_n,k,psi,Q,nodes,meshConfig,discretisationConsts);
    
    % Plot solutions at current time-step
%     solutionPlot(1) = subplot(2,3,1);
%     contourf(xNodes,zNodes,h_solved(:,:,end));
%     view(2)
%     colormap(solutionPlot(1),flipud(autumn))
%     shading interp;
%     colorbar
%     % water content
%     solutionPlot(2) = subplot(2,3,2);
%     contourf(xNodes,zNodes,psi_solved(:,:,end));
%     view(2)
%     colormap(solutionPlot(2),flipud(winter))
%     shading interp;
%     colorbar
%     % saturation
%     solutionPlot(3) = subplot(2,3,3);
%     contourf(xNodes,zNodes,S_solved(:,:,end));
%     view(2)
%     colormap(solutionPlot(3),cool)
%     shading interp;
%     colorbar
    % average water content (moisture)
    avgSatsMeasured(end+1) = sum(psi_h_n .* Deltas.xz)/(L1*L2);
%     solutionPlot(4) = subplot(2,3,[4,5,6]);
%     plot(t_vals,avgSatsMeasured,'b');
%     title("average water content at time " + num2str(t_vals(end)) + " (dt = " + num2str(dt) + ")");
%     drawnow;
end

% water table won't necessarily sit on the nodes - likely will be somewhere
% in the node control volume and need to determine it from pressure values


%% Show stuff to see where sim got up to
figure
psi_sat = psi(zeros(Nx*Nz,1))';
Psi_sat = sum((psi_sat(quadMats_p).*DVs.DV),2) ./ Deltas.xz;
Psi_sat_avg = sum(Psi_sat .* Deltas.xz)/(L1*L2);
avgSatsMeasured(end)/Psi_sat_avg

solutionPlot(1) = subplot(2,3,1);
contourf(xNodes,zNodes,h_solved(:,:,end));
view(2)
colormap(solutionPlot(1),flipud(autumn))
shading interp;
colorbar
title('Pressure Heads','Interpreter','Latex')
% water content
solutionPlot(2) = subplot(2,3,2);
contourf(xNodes,zNodes,psi_solved(:,:,end));
view(2)
colormap(solutionPlot(2),flipud(winter))
shading interp;
colorbar
title('Water Content','Interpreter','Latex')
% saturation
solutionPlot(3) = subplot(2,3,3);
contourf(xNodes,zNodes,S_solved(:,:,end));
view(2)
colormap(solutionPlot(3),flipud(parula))
shading interp;
colorbar
title('Saturation','Interpreter','Latex')
% average water content (moisture)
solutionPlot(4) = subplot(2,3,[4,5,6]);
plot(t_vals,avgSatsMeasured,'b');
ylim([0,0.5])
title('Average Water Content','Interpreter','Latex');
xlabel('Time (days)','Interpreter','Latex');
ylabel('$\psi_(avg)$','Interpreter','Latex')
drawnow;

%% Save solution to workspace
solution.t = t_vals;
solution.avgSats = avgSatsMeasured;
solution.h = h_solved;
solution.psi = psi_solved;
solution.S = S_solved;
solution.outflows = outflows;
solution.constants = discretisationConsts;
save('bs_rc_s90.mat', 'solution');

%% Load
load('bs_rc_s90.mat')
t_vals = solution.t;
avgSatsMeasured = solution.avgSats;
h_solved = solution.h;
psi_solved = solution.psi;
S_solved = solution.S;
outflows = solution.outflows;
discretisationConsts = solution.constants;


