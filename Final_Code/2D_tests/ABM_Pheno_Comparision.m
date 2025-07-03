clc;
close all;
clear all;

% profile on;

ModelParams = struct();

domainBoundary = struct();
domainBoundary.x_max =  50;
domainBoundary.y_max =  50;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 100;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 16;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 3;

% T-cells
ModelParams.P_A = 1.0; % Antigen gain rate
ModelParams.P_D = 0.5;  % Antigen loss rate
ModelParams.activatedAge  = 50;

dx = 0.5;
dy = dx;
dt = 0.1;
da = 1;

ModelParams.dx_ABM = dx;
ModelParams.dy_ABM = dy;
ModelParams.dt_ABM = dt;            % Time step
ModelParams.da_ABM = da;            % Time step
ModelParams.NumbAgents = 10000;

ModelParams.dx_PDE = dx;
ModelParams.dy_PDE = dy;
ModelParams.dt_PDE = dt;
ModelParams.da_PDE = da;

ModelParams.plot_traj = true;
recordRate = 20;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;


rndSeed = 2;

[u, ~, dCdx, dCdy, Ic, A_Ic, Inds_Ic_st, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams);
activation_proportion_pheno = zeros(params.nt,1);
x = linspace(0, params.Lx, params.Nx);
y = linspace(0, params.Ly, params.Ny);
% a = linspace(0, 1, params.Na);
a = linspace(0, params.La, params.Na);


[walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
activatedCells_ABM = zeros(params_ABM.nt,1);

% [walker_positions, walker_activation, U, C, Ap, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);

Ic = DC_model.Ic;
figure; 
subplot(1,2,1); surf(Ic, 'EdgeColor','none');
view(2)
subplot(1,2,2); surf(C,  'EdgeColor','none');
view(2)


%%

tic;

figure;
for n=0:params_ABM.nt
    n/params_ABM.nt;
     u = computePhenotypeModel(u, dCdx, dCdy, A, Inds_Ic_st, A_Ic, params);
    % activation_proportion_pheno(n) = squeeze(sum(u,[1,2]))'*a';
    activation_proportion_pheno(n+1) = squeeze(sum(u,[1,2]))'*a'/params.activatedAge;
    
    [walker_positions, walker_activation, ~, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
    activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*params_ABM.num_walkers);

    
    if false && mod(n,recordRate)==0
        clf;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ax1 = subplot(2,3,1);
        u_slice = squeeze(sum(u,3));
        surf(x, y, u_slice', 'EdgeColor','none');
        colorbar;
        view(2)
        ax1.XTick = [];
        ax1.YTick = [];
        xlabel('x', 'FontSize',14, 'Interpreter','latex');
        ylabel('y', 'FontSize',14, 'Interpreter','latex');
        title('PS-PDE T cell density', 'FontSize',14, 'Interpreter','latex');
        axis([0 params.Lx 0 params.Ly])
    
        av_a = zeros(params.Nx, params.Ny);
        for ix=1:params.Nx
            for iy=1:params.Ny
                av_a(ix,iy) = a*squeeze(u(ix,iy,:));
            end
        end
        av_a = av_a/(params.activatedAge*sum(u(:)));
        % av_a = av_a/(sum(u(:)));

    
        ax2 = subplot(2,3,2);
        surf(x, y, av_a', 'EdgeColor','none');
        colorbar;
        view(2)
        ax2.XTick = [];
        ax2.YTick = [];
        axis([0 params.Lx 0 params.Ly])
        xlabel('x', 'FontSize',14, 'Interpreter','latex');
        ylabel('y', 'FontSize',14, 'Interpreter','latex');
        title('PS-PDE antigen density', 'FontSize',14, 'Interpreter','latex');        
        subplot(2,3,3)
        plot(ModelParams.dt_PDE*(0:n),activation_proportion_pheno(1:(n+1)), 'LineWidth',2)
        axis([0 params.nt 0 1])

        subplot(2,3,6)
        % plot(params.activatedAge*a, 1/(params.activatedAge*params.da)*squeeze(sum(u,[1,2])), 'LineWidth',2)
        plot(a, 1/(params.da)*squeeze(sum(u,[1,2])), 'LineWidth',2)
        title('Antigen distribution')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        

        A_abm = zeros(params_ABM.Ny, params_ABM.Nx);
        U_abm = zeros(params_ABM.Ny, params_ABM.Nx);

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
            
        Ic = DC_model.Ic;
        U_plot = (U_abm/params_ABM.num_walkers);
        U_plot_1 = U_plot;
        % U_plot_1(Ic == 1) = nan;
        ax4 = subplot(2,3,4);
        hold on;
        surf(ax4,x, y, U_plot_1,'EdgeAlpha',0);
        hold off;
        xlabel('x', 'FontSize',14, 'Interpreter','latex');
        ylabel('y', 'FontSize',14, 'Interpreter','latex');
        title('ABM T cell density', 'FontSize',14, 'Interpreter','latex');
        colorbar;
        ax4.XTick = [];
        ax4.YTick = [];
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        view(0,90)

        tmp_act = (A_plot).*U_plot;
        % tmp_act(Ic == 1) = nan;
        ax5 = subplot(2,3,5);
        surf(ax5,x, y, tmp_act,'EdgeAlpha',0);
        xlabel('x', 'FontSize',14, 'Interpreter','latex');
        ylabel('y', 'FontSize',14, 'Interpreter','latex'); 
        title('ABM antigen density', 'FontSize',14, 'Interpreter','latex');        
        colorbar; % caxis([0 1]);
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        view(0,90)
        ax5.XTick = [];
        ax5.YTick = [];
        
        
        ax3 = subplot(2,3,3);
        hold on
        plot(ax3, ModelParams.dt_ABM*(0:n), activatedCells_ABM(1:(n+1)),'--','LineWidth',2.5);
        xlabel('Time', 'FontSize',14, 'Interpreter','latex'); 
        % ylabel('Total activation');
        axis([0 ModelParams.T_final 0 1])
        title('Proportion of activation', 'FontSize',14, 'Interpreter','latex')
    
        ax6 = subplot(2,3,6);
        hold on;
        [N_A, A_edges] = histcounts(walker_activation,-0.5:1:(params_ABM.activatedAge+0.5));
        plot(ax6, 0:1:params_ABM.activatedAge, N_A/ModelParams.NumbAgents,'--','LineWidth',2)
        title('Antigen distribution', 'FontSize',14, 'Interpreter','latex')


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




        
        sgtitle(['time = ', num2str(params.dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end


end

toc;

% profile viewer
% profile off;

