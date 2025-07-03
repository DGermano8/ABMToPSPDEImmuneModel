clear all;
% close all;

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];

addpath(genpath('UserFunctions'))

%%
numb_rands = 4;

ModelParams = struct();
rndSeed = 107;

domainBoundary = struct();
domainBoundary.x_max =  60;
domainBoundary.y_max =  60;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 1000;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 16;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 0.9; % Antigen gain rate
ModelParams.P_D = 0.9;  % Antigen loss rate
ModelParams.activatedAge  = 50;

ModelParams.dx_ABM = 1/3;
ModelParams.dy_ABM = 1/3;
ModelParams.dt_ABM = 0.025;        % Time step
ModelParams.NumbAgents = 1000;

ModelParams.dx_PDE = 1/3;
ModelParams.dy_PDE = 1/3;
ModelParams.dt_PDE = 0.025;

ModelParams.plot_traj = true;
recordRate = 25;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;

if (ModelParams.p_move*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    ModelParams.C_sens*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    (4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) == 1)

    disp('Choose smaller dx and dt for ABM')
else
    disp('Nice params!')
end


h = waitbar(0, 'Processing...'); % Create progress bar

%%

Cluster_Sizes = [1 2 3 4 8 16];

P_A_array = [0.1 0.5 0.9 0.1 0.5 0.9 0.1 0.5 0.9];
P_D_array = [0.1 0.1 0.1 0.5 0.5 0.5 0.9 0.9 0.9];

Activation_Cluster_ABM = zeros(length(P_A_array),length(Cluster_Sizes), numb_rands);
Activation_Cluster_PDE = zeros(length(P_A_array),length(Cluster_Sizes), numb_rands);

FullActivated_Cluster_ABM = zeros(length(P_A_array),length(Cluster_Sizes), numb_rands);
FullActivated_Cluster_PDE = zeros(length(P_A_array),length(Cluster_Sizes), numb_rands);

iter=0;

poolobj = gcp("nocreate"); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(4);
end


%%

tic;
for ips = 1:length(P_A_array)
    ModelParams.P_A = P_A_array(ips);
    ModelParams.P_D = P_D_array(ips);


    for cluster_i = 1:length(Cluster_Sizes)
    
        ModelParams.numberOfClusters = Cluster_Sizes(cluster_i);
        
        parfor jj=1:numb_rands
    
            [Cluster_Sizes(cluster_i) jj]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rndSeed = jj;
            % [U, A, U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params_PDE] = PDE_SetUp(rndSeed,domainBoundary,ModelParams);
            % activatedCells_PDE = zeros(params_PDE.nt+1,1);
            % for n = 0:params_PDE.nt
            %     [U, A, U_a, A_a] = computePDEModel(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy, params_PDE);
            %     tmp = (A+A_a);
            %     activatedCells_PDE(n+1) = ( sum( tmp(:) ) );
            % end
            % Activation_Cluster(cluster_i,jj) = activatedCells_PDE(end);
            % 
    
            [walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
            activatedCells_ABM = zeros(params_ABM.nt+1,1);
            activatedCells_T_ABM = zeros(params_ABM.nt+1,1);
    
            for n = 0:params_ABM.nt
                [walker_positions, walker_activation, U, C, Ap, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
                activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*sum(U(:)));
                activatedCells_T_ABM(n+1) = nnz(walker_activation==params_ABM.activatedAge)/length(walker_activation);
    
            end
            Activation_Cluster_ABM(ips,cluster_i,jj) = activatedCells_ABM(end);
            FullActivated_Cluster_ABM(ips,cluster_i,jj) = activatedCells_T_ABM(end);
            % plotABM(walker_positions, walker_activation, U, C, Ap, DCLingerTime, activatedCells_ABM, activatedCells_T_ABM, params_ABM, ModelParams, domainBoundary, DC_model);
            % drawnow;
    
            % [U, A,  U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params_PDE] = PDE_SetUp(rndSeed,domainBoundary,ModelParams);
            % activatedCells_PDE = zeros(params_PDE.nt+1,1);
            % activatedCells_T_PDE = zeros(params_PDE.nt+1,1);
            % for n = 0:params_PDE.nt
            %     [U, A, U_a, A_a] = computePDEModel(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy, params_PDE);
            %     tmp = (A+A_a);
            %     % tmp = (A+A_a)/params_PDE.activatedAge;
            %     % tmp = (A_a)/params_PDE.activatedAge;
            %     activatedCells_PDE(n+1) = ( sum( tmp(:) ) );
            %     activatedCells_T_PDE(n+1) = sum(A_a(:));
            % end
            % Activation_Cluster_PDE(ips,cluster_i,jj) = activatedCells_PDE(end);
            % FullActivated_Cluster_PDE(ips,cluster_i,jj) = activatedCells_T_PDE(end);
            % plotPDE(U, A, U_a, A_a, activatedCells_PDE, activatedCells_T_PDE, params_PDE, ModelParams, domainBoundary);
            % drawnow;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        end
        iter = iter + numb_rands;
        prog = (iter)/(numb_rands*length(Cluster_Sizes)*length(P_A_array));
    
        elapsedTime = toc;
        estimatedTotalTime = elapsedTime / prog;
        remainingTime = estimatedTotalTime - elapsedTime;
    
        if remainingTime < 120
            waitbar(prog, h, sprintf('Progress: %1.2f%%, Estimated time remaining: %1.1f seconds', 100*prog,  max(remainingTime, 0) )); % Update progress
        else
            waitbar(prog, h, sprintf('Progress: %1.2f%%, Estimated time remaining: %1.1f minutes', 100*prog,  max(remainingTime/60, 0) )); % Update progress
        end
    
    end
end
close(h);
total_time = toc;

%%
% 
% figure;
% subplot(2,2,1)
% violinplot(Activation_Cluster_ABM');
% xticklabels(["1"; "2"; "3"; "4"; "8"; "16"]')
% xlabel('mean # DC in cluster')
% ylabel('T-cell activation')
% title('ABM')
% axis([0 7 0 1])
% 
% subplot(2,2,2)
% violinplot(Activation_Cluster_PDE');
% xticklabels(["1"; "2"; "3"; "4"; "8"; "16"]')
% xlabel('mean # DC in cluster')
% ylabel('T-cell activation')
% title('PDE')
% axis([0 7 0 1])
% 
% 
% subplot(2,2,3)
% violinplot(FullActivated_Cluster_ABM');
% xticklabels(["1"; "2"; "3"; "4"; "8"; "16"]')
% xlabel('mean # DC in cluster')
% ylabel('T-cell activation')
% title('ABM')
% axis([0 7 0 1])
% 
% subplot(2,2,4)
% violinplot(FullActivated_Cluster_PDE');
% xticklabels(["1"; "2"; "3"; "4"; "8"; "16"]')
% xlabel('mean # DC in cluster')
% ylabel('T-cell activation')
% title('PDE')
% axis([0 7 0 1])
% 
% sgtitle(['\alpha+ = ', num2str(ModelParams.P_A), ', \alpha- = ', num2str(ModelParams.P_D), ', age = ', num2str(ModelParams.activatedAge)])

%%


figure;
for ips = 1:length(P_A_array)
    
    ModelParams.P_A = P_A_array(ips);
    ModelParams.P_D = P_D_array(ips);

    subplot(3,3,ips)
    ips_data = reshape(Activation_Cluster_ABM(ips,:,:), [length(Cluster_Sizes) numb_rands]);
    violinplot(ips_data');
    xticklabels(["1"; "2"; "3"; "4"; "8"; "16"]')
    % xlabel('mean # DC in cluster')
    % ylabel('T-cell activation')
    title(['\alpha+ = ', num2str(ModelParams.P_A), ', \alpha- = ', num2str(ModelParams.P_D)])
    axis([0 7 0 1])
    
end
sgtitle('Total Activation')

figure;
for ips = 1:length(P_A_array)
    
    ModelParams.P_A = P_A_array(ips);
    ModelParams.P_D = P_D_array(ips);

    ips_data = reshape(FullActivated_Cluster_ABM(ips,:,:), [length(Cluster_Sizes) numb_rands]);
    subplot(3,3,ips)
    violinplot(ips_data');
    xticklabels(["1"; "2"; "3"; "4"; "8"; "16"]')
    % xlabel('mean # DC in cluster')
    % ylabel('T-cell activation')
    title(['\alpha+ = ', num2str(ModelParams.P_A), ', \alpha- = ', num2str(ModelParams.P_D)])
    axis([0 7 0 1])
    
end
sgtitle('Full Activation')



%%

% figure;
% subplot(1,2,1)
% hold on;
% plot(mean(Activation_Cluster_ABM,2));
% plot(mean(Activation_Cluster_PDE,2));
% xlabel('mean # DC in cluster')
% ylabel('T-cell activation')
% 
% subplot(1,2,2)
% hold on;
% plot(mean(FullActivated_Cluster_ABM,2));
% plot(mean(FullActivated_Cluster_PDE,2));
% xlabel('mean # DC in cluster')
% ylabel('T-cell activation')

%% 

function plotABM(walker_positions, walker_activation, U, C, Ap, DCLingerTime, activatedCells_ABM, activatedCells_T_ABM, params_ABM, ModelParams, domainBoundary, DC_model)
    purple_colour = 1/255*[131, 96, 150];
    red_colour = 1/255*[237, 123, 123];
    yellow_colour = 1/255*[240, 184, 110];
    blue_colour = 1/255*[106 123 237];
    green_colour = 1/255*[77 150 138];
    back_colour = 1/255*[56 56 56];
    
    cmap = myColour3Gradient(255, yellow_colour, red_colour, purple_colour);
    


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

    figure;
    x = linspace(0, params_ABM.Lx, params_ABM.Nx);
    y = linspace(0, params_ABM.Ly, params_ABM.Ny);
    [X, Y] = meshgrid(x, y);
    total_density = sum(U(:));
    title_string = ['Total activation = ', num2str(100*activatedCells_ABM(params_ABM.nt),'%.1f')];
    % sgtitle(title_string,'Color',1/255*[217 217 217],'FontSize',18);
    
    Ic = DC_model.Ic;
    U_plot = (U_abm/params_ABM.num_walkers);
    U_plot_1 = U_plot;
    U_plot_1(Ic == 1) = nan;
    ax1 = subplot(2,3,1);
    hold on;
    surf(ax1,X, Y, U_plot_1,'EdgeAlpha',0);
    hold off;
    xlabel('x'); ylabel('y'); title('U');
    c1 = colorbar;
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    % colormap(ax1,cmap2);
    caxis([0 max(U_plot(:))])
    view(0,90)
    ax1.XTick = [0:10:50];
    ax1.YTick = [0:10:50];
    c1.Label.String = 'T-cell density';
    c1.Label.FontSize = 12;
    colormap(cmap);
    
    tmp_act = (A_plot).*U_plot;
    tmp_act(Ic == 1) = nan;
    ax2 = subplot(2,3,2);
    surf(ax2,X, Y, tmp_act,'EdgeAlpha',0);
    xlabel('x'); ylabel('y'); title('A');
    c2 = colorbar; % caxis([0 1]);
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    view(0,90)
    ax2.XTick = [0:10:50];
    ax2.YTick = [0:10:50];
    c2.Label.String = 'T-cell activation';
    c2.Label.FontSize = 12;
    
    subplot(2,3,3)
    histogram(walker_activation)
    title(title_string,'FontSize',12);
    
    subplot(2,3,6);
    hold on
    plot( 0:ModelParams.dt_ABM:ModelParams.T_final, activatedCells_ABM,'LineWidth',2.5);
    plot( 0:ModelParams.dt_ABM:ModelParams.T_final, activatedCells_T_ABM,'LineWidth',2.5);
    xlabel('Time'); ylabel('Total activation');
    axis([0 ModelParams.T_final 0 1])
    drawnow;


end

function plotPDE(U, A, U_a, A_a, activatedCells_PDE, activatedCells_T_PDE, params_PDE, ModelParams, domainBoundary)

    purple_colour = 1/255*[131, 96, 150];
    red_colour = 1/255*[237, 123, 123];
    yellow_colour = 1/255*[240, 184, 110];
    blue_colour = 1/255*[106 123 237];
    green_colour = 1/255*[77 150 138];
    back_colour = 1/255*[56 56 56];
    
    cmap = myColour3Gradient(255, yellow_colour, red_colour, purple_colour);
   
    x = linspace(0, params_PDE.Lx, params_PDE.Nx);
    y = linspace(0, params_PDE.Ly, params_PDE.Ny);
    [X, Y] = meshgrid(x, y);
    
    total_density = sum(U(:));
    title_string = ['Total activation = ', num2str(100*activatedCells_PDE(end),'%.1f')];
    
    subplot(2,3,4)
    surf(X, Y, (U/(sum(U(:)))),'EdgeAlpha',0);
    xlabel('x'); ylabel('y'); title('U');
    colorbar;
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    view(0,90)
    
    subplot(2,3,5)
    % surf(X, Y, (A)/(ModelParams.activatedAge),'EdgeAlpha',0);
    surf(X, Y, (A+A_a),'EdgeAlpha',0);
    xlabel('x'); ylabel('y'); title('A');
    colorbar; 
    axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
    view(0,90)
    
    subplot(2,3,6);
    hold on;
    plot(0:ModelParams.dt_PDE:ModelParams.T_final, activatedCells_PDE, '--' ,'LineWidth',2.5)
    plot(0:ModelParams.dt_PDE:ModelParams.T_final, activatedCells_T_PDE, ':' ,'LineWidth',2.5)
    axis([0 ModelParams.T_final 0 1])
    title(title_string,'FontSize',12);
    drawnow;
    hold off;
end


