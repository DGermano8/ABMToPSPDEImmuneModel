% clc;
close all;
clear all;


purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];

%%
% close all;
% profile on;

ModelParams = struct();

domainBoundary = struct();
domainBoundary.x_max =  40;
domainBoundary.y_max =  40;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 600;
ModelParams.p_move = 0.2;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 7;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 0.8; % Antigen gain rate
ModelParams.P_D = 0.4;  % Antigen loss rate
ModelParams.activatedAge  = 50;

dx = 0.25;
dy = dx;
dt = 0.005;
da = 0.25;

ModelParams.dx_ABM = dx;
ModelParams.dy_ABM = dy;
ModelParams.dt_ABM = dt;            % Time step
ModelParams.da_ABM = da;            % Time step
ModelParams.NumbAgents = 1;

ModelParams.dx_PDE = dx;
ModelParams.dy_PDE = dy;
ModelParams.dt_PDE = dt;
ModelParams.da_PDE = da;

ModelParams.plot_traj = true;
recordRate = 50/dt;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;
step_plot = 5;

% rndSeed = 6; % for 1 DC and domain -> 9
rndSeed = 7;% for 1 DC and domain -> 9

[(ModelParams.p_move)*(4*dt)/(dx.^2) (ModelParams.C_sens)*(4*dt)/(dx.^2) (ModelParams.P_A/da)*dt (ModelParams.P_A/ModelParams.activatedAge)*dt (ModelParams.P_D/da)*dt ]



[walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
activatedCells_ABM = zeros(params_ABM.nt,1);
% seed = rndSeed;

x = linspace(0, params_ABM.Lx, params_ABM.Nx);
y = linspace(0, params_ABM.Ly, params_ABM.Ny);
a = linspace(0, params_ABM.La, params_ABM.Na);

tmp_Ic = DC_model.Ic;
tmp_Ic(tmp_Ic == 0) = nan;

%%
figure;

surf(x,y,tmp_Ic, 'FaceColor','black','EdgeColor','none')
axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
xlabel('x', 'FontSize',12,'Interpreter','latex')
ylabel('y', 'FontSize',12,'Interpreter','latex')
xticks([0:10:40])
yticks([0:10:40])
view(2);