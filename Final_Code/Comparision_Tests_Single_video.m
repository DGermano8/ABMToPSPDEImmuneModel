clear all;
% close all;

addpath(genpath('UserFunctions'))



RecordVideo = true;
recordRate = 10;

f = figure;
f.Position = [106 360 960 580];

if RecordVideo
    fps = 15;
    fileName = 'test6.mp4';
    
    vidfile = VideoWriter(fileName,'MPEG-4');
    vidfile.FrameRate = fps;
    open(vidfile);
    numberOfSteps = 0;
end

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];
white_colour = 1/255*[217 217 217];

cmap = myColour2Gradient(255, red_colour, blue_colour);

plot_abm = true;
plot_pde = true;

%%

ModelParams = struct();
rndSeed = 106;

domainBoundary = struct();
domainBoundary.x_max =  7;
domainBoundary.y_max =  7;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 1000;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 1;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 1;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 0.2; % Antigen gain rate
ModelParams.P_D = 0.2;  % Antigen loss rate
ModelParams.activatedAge  = 50;

% ModelParams.dx_ABM = 0.25;
% ModelParams.dy_ABM = 0.25;
% ModelParams.dt_ABM = 0.01;            % Time step
% ModelParams.NumbAgents = 5000;
% 
% ModelParams.dx_PDE = 0.25;
% ModelParams.dy_PDE = 0.25;
% ModelParams.dt_PDE = 0.01;

ModelParams.dx_ABM = 1/3;
ModelParams.dy_ABM = 1/3;
ModelParams.dt_ABM = 0.025;   % Time step
ModelParams.NumbAgents = 10000;

ModelParams.dx_PDE = 1/3;
ModelParams.dy_PDE = 1/3;
ModelParams.dt_PDE = 0.025;

ModelParams.plot_traj = true;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;

if (ModelParams.p_move*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    ModelParams.C_sens*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    (4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) > 1)
    clc;
    disp('Choose smaller dx and dt for ABM')
else
    clc;
    disp('Nice params!')
    disp([' ps = ', num2str((4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2))] )

end
% h = waitbar(0, 'Processing...'); % Create progress bar

%%
tic;
[walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
activatedCells_ABM = zeros(params_ABM.nt+1,1);
activatedCells_T_ABM = zeros(params_ABM.nt+1,1);


[U, A, U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params_PDE] = PDE_SetUp(rndSeed,domainBoundary,ModelParams);
activatedCells_PDE = zeros(params_PDE.nt+1,1);
activatedCells_T_PDE = zeros(params_PDE.nt+1,1);


for n = 0:params_ABM.nt
    
    [walker_positions, walker_activation, U_abm, C_abm, Ap_abm, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
    % [walker_positions, walker_activation, U_abm, C_abm, Ap_abm, DCLingerTime] = computeABMModel_Vectorized_fixedProb(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
    activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*length(walker_activation));
    activatedCells_T_ABM(n+1) = nnz(walker_activation==params_ABM.activatedAge)/length(walker_activation);

    
    [U, A, U_a, A_a] = computePDEModel(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy, params_PDE);
    tmp = (A+A_a);
    activatedCells_PDE(n+1) = ( sum( tmp(:) ) )/(ModelParams.activatedAge);
    activatedCells_T_PDE(n+1) = ( sum( A_a(:) ) )/(ModelParams.activatedAge);
    
    if mod(n, (recordRate/ModelParams.dt_PDE) ) == 0

        clf;
        hold on;

        Ic = DC_model.Ic;
        U_p = U(Ic == 0);
        U_caxis = [min(U_p(:)) max(U_p(:))];
        A_T = (A+A_a)/(ModelParams.activatedAge);
        A_T = A_T(Ic == 0);
        A_caxis = [min(A_T(:)) max(max(A_T(:)), 10^-6)];

        plotABM(n, walker_positions, walker_activation, U_abm, C_abm, Ap_abm, DCLingerTime, activatedCells_ABM, activatedCells_T_ABM, params_ABM, ModelParams, domainBoundary, DC_model, U_caxis, A_caxis);

        plotPDE(n, U, A/(ModelParams.activatedAge), U_a, A_a/(ModelParams.activatedAge), activatedCells_PDE, activatedCells_T_PDE, params_PDE, ModelParams, domainBoundary, Inds_Ic.Ic_vals, U_caxis, A_caxis);

        sgtitle(['Time = ' num2str(n*ModelParams.dt_PDE)], 'Color',white_colour, 'FontSize',18, 'Interpreter','latex');
        
        histogram(walker_activation,0:ModelParams.activatedAge+1)
        
        MakeDark();
        drawnow;

        if RecordVideo
            numberOfSteps = numberOfSteps + 1;
            F(numberOfSteps) = getframe(gcf);
            writeVideo(vidfile,F(numberOfSteps));
        end
    end
end

if RecordVideo
    close(vidfile)
end


%%
Ic = DC_model.Ic;
U_p = U(Ic == 0);
U_caxis = [min(U_p(:)) max(U_p(:))];
A_T = (A+A_a)/(ModelParams.activatedAge);
A_T = A_T(Ic == 0);
A_caxis = [min(A_T(:)) max(max(A_T(:)), 10^-6)];

plotABM(n, walker_positions, walker_activation, U_abm, C_abm, Ap_abm, DCLingerTime, activatedCells_ABM, activatedCells_T_ABM, params_ABM, ModelParams, domainBoundary, DC_model, U_caxis, A_caxis);

plotPDE(n, U, A/(ModelParams.activatedAge), U_a, A_a/(ModelParams.activatedAge), activatedCells_PDE, activatedCells_T_PDE, params_PDE, ModelParams, domainBoundary, Inds_Ic.Ic_vals, U_caxis, A_caxis);

sgtitle(['Time = ' num2str(n*ModelParams.dt_PDE)], 'Color',white_colour, 'FontSize',18, 'Interpreter','latex');

histogram(walker_activation,0:ModelParams.activatedAge+1)

MakeDark();
drawnow;

%%


% figure;
% NumberDCsToActivate = zeros(ModelParams.NumDCs,1);
% for ii=1:size(DCLingerTime,1)
% 
%     if walker_activation(ii) == ModelParams.activatedAge
%         ii_NumberDCsToActivate = nnz(DCLingerTime(ii,:));
%         if ii_NumberDCsToActivate > 0
%             NumberDCsToActivate(ii_NumberDCsToActivate) = NumberDCsToActivate(ii_NumberDCsToActivate) + 1;
%         end
%     end
% end
% 
% % subplot(2,2,4)
% % Plot the histogram using bar
% if ModelParams.P_A < 0.5
%     bar(1:ModelParams.NumDCs, NumberDCsToActivate, 'FaceColor', red_colour);
% else
%     bar(1:ModelParams.NumDCs, NumberDCsToActivate, 'FaceColor', blue_colour);
% end
% % Add labels and title
% xlabel('Number of DCs', 'FontSize',12, 'Interpreter','latex');
% ylabel('Counts', 'FontSize',12, 'Interpreter','latex');
% axis([0 ModelParams.NumDCs 0 1.05*max(NumberDCsToActivate)])
% 
% MakeDark();

%%

function plotABM(n, walker_positions, walker_activation, U, C, Ap, DCLingerTime, activatedCells_ABM, activatedCells_T_ABM, params_ABM, ModelParams, domainBoundary, DC_model, U_caxis, A_caxis)
    purple_colour = 1/255*[131, 96, 150];
    red_colour = 1/255*[237, 123, 123];
    yellow_colour = 1/255*[240, 184, 110];
    blue_colour = 1/255*[106 123 237];
    green_colour = 1/255*[77 150 138];
    back_colour = 1/255*[56 56 56];
    
    cmap = myColour3Gradient(255, red_colour, purple_colour, blue_colour);
    


    A_abm = zeros(size(U));
    U_abm = zeros(size(U));
    for i = 1:params_ABM.num_walkers
        id_x_i = round(walker_positions(i,1)/params_ABM.dx) + 1;
        id_y_i = round(walker_positions(i,2)/params_ABM.dy) + 1;
    
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
    % total_density = sum(U(:));
    % title_string = ['Total activation = ', num2str(100*activatedCells_ABM(params_ABM.nt),'%.1f')];
    % sgtitle(title_string,'Color',1/255*[217 217 217],'FontSize',18);
    
    Ic = DC_model.Ic;
    U_plot = (U_abm/params_ABM.num_walkers);
    U_plot_1 = U_plot;
    U_plot_1(Ic == 1) = nan;
    ax1 = subplot(2,3,1);
    hold on;
    surf(ax1,X, Y, U_plot_1,'EdgeAlpha',0);
    hold off;
    xlabel('x', 'FontSize',14, 'Interpreter','latex');
    ylabel('y', 'FontSize',14, 'Interpreter','latex');
    title('T cell density', 'FontSize',14, 'Interpreter','latex');
    colorbar;
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    % colormap(ax1,cmap2);
    caxis(U_caxis)
    view(0,90)
    ax1.XTick = [];
    ax1.YTick = [];
    % ax1.XTick = [0 3.5 7];
    % ax1.YTick = [0 3.5 7];
    % c1.Label.String = 'T-cell density';
    % c1.Label.FontSize = 12;
    colormap(cmap);
    
    tmp_act = (A_plot).*U_plot;
    tmp_act(Ic == 1) = nan;
    ax2 = subplot(2,3,2);
    surf(ax2,X, Y, tmp_act,'EdgeAlpha',0);
    xlabel('x', 'FontSize',14, 'Interpreter','latex');
    ylabel('y', 'FontSize',14, 'Interpreter','latex'); 
    title('Antigen density', 'FontSize',14, 'Interpreter','latex');
    colorbar; % caxis([0 1]);
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    view(0,90)
    ax2.XTick = [];
    ax2.YTick = [];
    % ax2.XTick = [0 3.5 7];
    % ax2.YTick = [0 3.5 7];
    caxis(A_caxis)
    % c2.Label.String = 'T-cell activation';
    % c2.Label.FontSize = 12;
    
    % subplot(2,3,3)
    % histogram(walker_activation)
    % title(title_string,'FontSize',12);
    
    subplot(2,3,3);
    hold on
    plot( ModelParams.dt_ABM*(0:n), activatedCells_ABM(1:(n+1)),'LineWidth',2.5);
    xlabel('Time', 'FontSize',14, 'Interpreter','latex'); 
    % ylabel('Total activation');
    axis([0 ModelParams.T_final 0 1])
    title('Proportion of activation', 'FontSize',14, 'Interpreter','latex')

    subplot(2,3,6)
    plot( ModelParams.dt_ABM*(0:n), activatedCells_T_ABM(1:(n+1)),'LineWidth',2.5);
    xlabel('Time', 'FontSize',14, 'Interpreter','latex');
    % label('Total activation');
    axis([0 ModelParams.T_final 0 1])
    title('Proportion of activated T cells', 'FontSize',14, 'Interpreter','latex')

end

function plotPDE(n, U, A, U_a, A_a, activatedCells_PDE, activatedCells_T_PDE, params_PDE, ModelParams, domainBoundary, Inds_Ic, U_caxis, A_caxis)

    purple_colour = 1/255*[131, 96, 150];
    red_colour = 1/255*[237, 123, 123];
    yellow_colour = 1/255*[240, 184, 110];
    blue_colour = 1/255*[106 123 237];
    green_colour = 1/255*[77 150 138];
    back_colour = 1/255*[56 56 56];
    
    cmap = myColour3Gradient(255, red_colour, purple_colour, blue_colour);
   
    x = linspace(0, params_PDE.Lx, params_PDE.Nx);
    y = linspace(0, params_PDE.Ly, params_PDE.Ny);
    [X, Y] = meshgrid(x, y);
    
    % total_density = sum(U(:));
    % title_string = ['Total activation = ', num2str(100*activatedCells_PDE(end),'%.1f')];
    
    ax3 = subplot(2,3,4);
    U(Inds_Ic) = nan;

    surf(ax3, X, Y, U,'EdgeAlpha',0);
    xlabel('x', 'FontSize',14, 'Interpreter','latex');
    ylabel('y', 'FontSize',14, 'Interpreter','latex');
    title('T cell density', 'FontSize',14, 'Interpreter','latex');
    colorbar;
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    view(0,90)
    ax3.XTick = [];
    ax3.YTick = [];
    % ax3.XTick = [0 3.5 7];
    % ax3.YTick = [0 3.5 7];
    caxis(U_caxis)
    colormap(cmap);

    ax4 = subplot(2,3,5);
    % surf(X, Y, (A)/(ModelParams.activatedAge),'EdgeAlpha',0);
    A(Inds_Ic) = nan;
    A_a(Inds_Ic) = nan;
    
    A_T = (A+A_a);
    surf(ax4, X, Y, A_T,'EdgeAlpha',0);
    xlabel('x', 'FontSize',14, 'Interpreter','latex');
    ylabel('y', 'FontSize',14, 'Interpreter','latex');
    title('Antigen density', 'FontSize',14, 'Interpreter','latex')
    colorbar;
    caxis(A_caxis)
    ax4.XTick = [];
    ax4.YTick = [];
    % ax4.XTick = [0 3.5 7];
    % ax4.YTick = [0 3.5 7];
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    view(0,90)
    
    subplot(2,3,3);
    hold on;
    plot(ModelParams.dt_PDE*(0:n), activatedCells_PDE(1:(n+1)), '--' ,'LineWidth',2.5)
    axis([0 ModelParams.T_final 0 1])
    % title(title_string,'FontSize',12);
    title('Proportion of activation', 'FontSize',14, 'Interpreter','latex')

    subplot(2,3,6);
    hold on;
    plot(ModelParams.dt_PDE*(0:n), activatedCells_T_PDE(1:(n+1)), '--' ,'LineWidth',2.5)
    axis([0 ModelParams.T_final 0 1])
    % title(title_string,'FontSize',12);
    title('Proportion of activated T cells', 'FontSize',14, 'Interpreter','latex')
    hold off;
end


