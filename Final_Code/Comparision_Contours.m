% clear all;
% close all;
addpath(genpath('UserFunctions'))
figure;
 
 
purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56]; 
 
%%

ModelParams = struct();
rndSeed = 108;

domainBoundary = struct();
domainBoundary.x_max =  7;
domainBoundary.y_max =  7;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 100;
ModelParams.p_move = 0.5;
 
% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;
 
% DCs
ModelParams.NumDCs = 1;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;
 
% T-cells
ModelParams.activatedAge  = 100;
ModelParams.P_A = 0.5; % Antigen gain rate
ModelParams.P_D = 0.5;  % Antigen loss rate

ModelParams.dx_ABM = 0.25;
ModelParams.dy_ABM = 0.25;
ModelParams.dt_ABM = 0.01;            % Time step
ModelParams.NumbAgents = 1000000;
 
ModelParams.dx_PDE = 0.25;
ModelParams.dy_PDE = 0.25;
ModelParams.dt_PDE = 0.01;

ModelParams.plot_traj = true;
recordRate = 25;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;

plt_timer = 0;
plot_times = [0 5 10 25 50 100]/ModelParams.dt_ABM;
cmap = myColour3Gradient(100+1, yellow_colour, red_colour, purple_colour);
cmap = myColour3Gradient(6, yellow_colour, red_colour, purple_colour);
cmap = custom_colormap(yellow_colour, red_colour, purple_colour, 6)

 
if (ModelParams.p_move*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    ModelParams.C_sens*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    (4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) > 1)
    clc;
    disp('Choose smaller dx and dt for ABM')
else
    clc;
    disp('Nice params!')
 
end
disp([' ps = ', num2str((4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2))] )
h = waitbar(0, 'Processing...'); % Create progress bar
 
%%
tic;
[walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
activatedCells_ABM = zeros(params_ABM.nt+1,1);
for n = 0:params_ABM.nt
 
    [walker_positions, walker_activation, U, C, Ap, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
    activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*sum(U(:)));
    % activatedCells_ABM(n+1) = sum(( walker_activation == ModelParams.activatedAge))/(length(walker_activation));
    if mod(n, (1/ModelParams.dt_ABM) ) == 0
    progress = n/params_ABM.nt;
    elapsedTime = toc;
    estimatedTotalTime = elapsedTime / progress;
    remainingTime = estimatedTotalTime - elapsedTime;
 
    waitbar(n/params_ABM.nt, h, sprintf('Progress: %1.2f%%, Estimated time remaining: %1.1f seconds', 100*progress,  max(remainingTime, 0) )); % Update progress
    end

    if ismember(n,plot_times)
        plt_timer = plt_timer + 1;
        A_abm = zeros(size(U));
        U_abm = zeros(size(U));
        for i = 1:params_ABM.num_walkers
            id_x_i = walker_positions(i,1)/params_ABM.dx + 1;
            id_y_i = walker_positions(i,2)/params_ABM.dy + 1;
         
            U_abm(id_y_i,id_x_i) = U_abm(id_y_i,id_x_i) + 1;
            A_abm(id_y_i,id_x_i) = A_abm(id_y_i,id_x_i) + walker_activation(i);
        end
        A_plot = zeros(size(U_abm));
        for ii=1:size(U_abm,1)
            for jj=1:size(U_abm,2)
                if U_abm(jj,ii) > 0 
                    A_plot(jj,ii) = A_abm(jj,ii)/(params_ABM.activatedAge*U_abm(jj,ii));
                end
            end
        end
         
        x = linspace(0, params_ABM.Lx, params_ABM.Nx);
        y = linspace(0, params_ABM.Ly, params_ABM.Ny);
        [X, Y] = meshgrid(x, y);
        total_density = sum(U(:));
        title_string = ['Total activation = ', num2str(100*activatedCells_ABM(params_ABM.nt),'%.1f')];
        % sgtitle(title_string,'Color',1/255*[217 217 217],'FontSize',18);
         
        U_plot = (U_abm/params_ABM.num_walkers);
        subplot(1,3,1);
        hold on;
        % surf(X, Y, U_plot,'EdgeAlpha',0);
        % contour(X, Y, U_plot,[1.5*10^(-3) 1.5*10^(-3)], 'color', cmap(1+n*ModelParams.dt_ABM,:));
        contour(X, Y, U_plot,[1.5*10^(-3) 1.5*10^(-3)], 'color', cmap(plt_timer,:), 'linewidth',4);
        hold off;
        xlabel('x', 'Interpreter','latex', 'FontSize',14)
        ylabel('y', 'Interpreter','latex', 'FontSize',14)
        title('T cell density', 'Interpreter','latex', 'FontSize',16)
        % c1 = colorbar;
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        % colormap(ax1,cmap2);
        caxis([0 max(U_plot(:))])
        view(0,90)
        XTick = [0:10:50];
        YTick = [0:10:50];
        % c1.Label.String = 'T-cell density';
        % c1.Label.FontSize = 12;
        colormap(cmap);
         
        subplot(1,3,2);
        hold on;
        % surf(X, Y, (A_plot).*U_plot,'EdgeAlpha',0);
        % contour(X, Y, A_plot.*U_plot, [3*10^(-4) 3*10^(-4)], 'color', cmap(1+n*ModelParams.dt_ABM,:))
        contour(X, Y, A_plot.*U_plot, [0.15*10^(-4) 0.15*10^(-4)], 'color', cmap(plt_timer,:), 'linewidth',4)
        xlabel('x', 'Interpreter','latex', 'FontSize',14)
        ylabel('y', 'Interpreter','latex', 'FontSize',14)
        title('Antigen density', 'Interpreter','latex', 'FontSize',16)
        % c2 = colorbar; % caxis([0 1]);
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        view(0,90)
        XTick = [0:10:50];
        YTick = [0:10:50];
        % c2.Label.String = 'T-cell activation';
        % c2.Label.FontSize = 12;
         
    end
end
% close(h);


subplot(1,3,3)
hold on;
plot( 0:ModelParams.dt_ABM:ModelParams.T_final, activatedCells_ABM, ...
     '-', 'LineWidth',3, 'Color',blue_colour)
xlabel('Time'); ylabel('Total activation');
axis([0 ModelParams.T_final 0 1])


plt_timer = 0;
[U, A, U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params_PDE] = PDE_SetUp(rndSeed,domainBoundary,ModelParams);
activatedCells_PDE = zeros(params_PDE.nt+1,1);
for n = 0:params_PDE.nt
    if mod(n, (1/ModelParams.dt_PDE) ) == 0
        waitbar(n/params_PDE.nt, h, sprintf('Progress: %1.2f%%', 100*n/params_PDE.nt)); % Update progress
    end
    [U, A, U_a, A_a] = computePDEModel(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy, params_PDE);
    % tmp = (A)/params_PDE.activatedAge;
    % tmp = A.*U;
    tmp = (A+A_a);
    % tmp = (A_a)/params_PDE.activatedAge;
    % tmp = sum(U_a(:));
    activatedCells_PDE(n+1) = ( sum( tmp(:) ) );

    if ismember(n,plot_times)
        plt_timer = plt_timer + 1;
        x = linspace(0, params_PDE.Lx, params_PDE.Nx);
        y = linspace(0, params_PDE.Ly, params_PDE.Ny);
        [X, Y] = meshgrid(x, y);
        
        total_density = sum(U(:));
        title_string = ['Total activation = ', num2str(100*activatedCells_PDE(n+1),'%.1f')];
         
        subplot(1,3,1)
        hold on;
        % surf(X, Y, (U/(sum(U(:)))),'EdgeAlpha',0);
        % contour(X, Y, (U/(sum(U(:)))),[1.5*10^(-3) 1.5*10^(-3)], 'color', cmap(1+n*ModelParams.dt_PDE,:))
        % contour(X, Y, (U/(sum(U(:)))),[1.5*10^(-3) 1.5*10^(-3)], 'color', cmap(plt_timer,:), 'linewidth',2.5)
        contour(X, Y, (U/(sum(U(:)))),[1.5*10^(-3) 1.5*10^(-3)] ,  ':', 'color', [0 0 0], 'linewidth',2)
        xlabel('x', 'Interpreter','latex', 'FontSize',14)
        ylabel('y', 'Interpreter','latex', 'FontSize',14)
        title('T cell density', 'Interpreter','latex', 'FontSize',16)
        % colorbar;
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        view(0,90)
        
        subplot(1,3,2)
        hold on;
        % surf(X, Y, (A)/(ModelParams.activatedAge),'EdgeAlpha',0);
        % surf(X, Y, (A+A_a),'EdgeAlpha',0);
        % contour(X, Y, (A+A_a), [3*10^(-4) 3*10^(-4)], 'color', cmap(1+n*ModelParams.dt_PDE,:))
        % contour(X, Y, (A+A_a), [0.15*10^(-4) 0.15*10^(-4)], 'color', cmap(plt_timer,:), 'linewidth',2.5)
        contour(X, Y, (A+A_a), [0.15*10^(-4) 0.15*10^(-4)],  ':', 'color', [0 0 0], 'linewidth',2)
        xlabel('x', 'Interpreter','latex', 'FontSize',14)
        ylabel('y', 'Interpreter','latex', 'FontSize',14)
        title('Antigen density', 'Interpreter','latex', 'FontSize',16)
        % colorbar; 
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        view(0,90)        
    end
end
close(h);

subplot(1,3,1)
xticks([0 3.5 7])
yticks([0 3.5 7])

subplot(1,3,2)
xticks([0 3.5 7])
yticks([0 3.5 7])
 
subplot(1,3,3)
hold on;
plot(0:ModelParams.dt_PDE:ModelParams.T_final, activatedCells_PDE, '--', 'LineWidth',2.5, 'Color',red_colour)
axis([0 ModelParams.T_final 0 1])
xlabel('t', 'Interpreter','latex', 'FontSize',14)
ylabel('T cell activation', 'Interpreter','latex', 'FontSize',14)
title('Proportion of activated T cells', 'Interpreter','latex', 'FontSize',16)
legend('ABM','PDE', 'Interpreter','latex', 'FontSize',14)

xticks([0 50 100])
yticks([0 0.5 1])
drawnow;

f = gcf;
f.Color = [1 1 1];
% export_fig 1D_NoTaxis_NoLoss.png -m2.5


%%
 
% Define two different x-meshes
x1 = 0:ModelParams.dt_PDE:ModelParams.T_final'; % First x-mesh
y1 = activatedCells_PDE'; % First curve
 
x2 = 0:ModelParams.dt_ABM:ModelParams.T_final'; % Second x-mesh
y2 = activatedCells_ABM'; % Second curve
 
% Interpolate y2 onto x1's mesh
y2_interp = interp1(x2, y2, x1, 'previous'); % Linear interpolation
 
    
% Compute the absolute area between the curves
area_between = (trapz(x1, abs(y1 - y2_interp)));
area_between = (area_between)/ModelParams.T_final*100;
 
 
disp(['Area between curves: ', num2str(area_between)]);