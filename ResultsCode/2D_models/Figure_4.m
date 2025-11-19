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

walker_path = zeros(params_ABM.nt,2);
cmap = myColour3Gradient(params_ABM.Na, yellow_colour, red_colour, purple_colour); % Use parula colormap

tic;

figure;
for n=0:params_ABM.nt
    % seed = seed + 1;
    n/params_ABM.nt;
    
    
    % Call the C MEX function
    [walker_positions, walker_activation, C, DCLingerTime] = ...
        computeABMModel_mex( ...
            walker_positions, ...       % 1
            walker_activation, ...      % 2
            C, ...                      % 3
            DCLingerTime, ...           % 4
            DC_model.Ic, ...            % 5
            DC_model.BoundaryDC, ...    % 6
            params_ABM.num_walkers, ... % 7
            params_ABM.dx, ...          % 8
            params_ABM.dy, ...          % 9
            params_ABM.Nx, ...          % 10
            params_ABM.Ny, ...          % 11
            params_ABM.p_move, ...      % 12
            params_ABM.C_chi, ...       % 13
            params_ABM.activatedAge, ...% 14
            params_ABM.P_A, ...         % 15
            params_ABM.P_D, ...         % 16
            params_ABM.dt, ...          % 17
            params_ABM.da, ...          % 18
            ModelParams.NumDCs ...
        );
    walker_path(n+1,:) = walker_positions;
    activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*params_ABM.num_walkers);

    
    if true && mod(n,recordRate)==0
        clf;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        subplot(1,2,1)
        hold on;
        xx = walker_path(1:step_plot:(n+1),1)';
        yy = walker_path(1:step_plot:(n+1),2)';
        
        % Define a color value for each point (e.g., based on x)
        c = activatedCells_ABM(1:step_plot:(n+1))';
        
        % Plot line with gradient
        surface([xx; xx], [yy; yy], [zeros(size(xx)); zeros(size(xx))], ...
            [c; c], 'facecolor', 'none', 'edgecolor', 'interp', 'linewidth', 2.5);
        
        colormap(cmap)      % choose any colormap (e.g. 'parula', 'hot', 'cool')
        cb = colorbar;              % Create colorbar and store its handle
        cb.Label.String = 'Stimulation Level';
        cb.Label.FontSize = 12;
        cb.Label.Interpreter = 'latex';
        caxis([0,1])
        

        surf(x,y,tmp_Ic, 'FaceColor','black','EdgeColor','none')
        scatter3(walker_positions(1),walker_positions(2),0,75,'filled','k','square')
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        xlabel('x', 'FontSize',12,'Interpreter','latex')
        ylabel('y', 'FontSize',12,'Interpreter','latex')
        title('ABM T cell - Dendritic cell interaction', 'FontSize',14,'Interpreter','latex')
        xticks([0:10:40])
        yticks([0:10:40])

        subplot(1,2,2)
        % plot(dt*(1:(n+1)),activatedCells_ABM(1:(n+1)))
        xxx = dt*(1:step_plot:(n+1));
        yyy = activatedCells_ABM(1:step_plot:(n+1))';
        
        % Define a color value for each point (e.g., based on x)
        ccc = activatedCells_ABM(1:step_plot:(n+1))';
        
        % Plot line with gradient
        surface([xxx; xxx], [yyy; yyy], [zeros(size(xxx)); zeros(size(xxx))], ...
            [ccc; ccc], 'facecolor', 'none', 'edgecolor', 'interp', 'linewidth', 3);

        axis([0 dt*params_ABM.nt 0 1])
        xlabel('Time (hrs)', 'FontSize',12,'Interpreter','latex')
        ylabel('T cell stimulation level', 'FontSize',12,'Interpreter','latex')
        title('ABM T cell stimulation level', 'FontSize',14,'Interpreter','latex')
        yticks([0:0.2:1])
        % sgtitle(['time = ', num2str(params_ABM.dt*n)])
        drawnow;
        
    end

end

toc;



