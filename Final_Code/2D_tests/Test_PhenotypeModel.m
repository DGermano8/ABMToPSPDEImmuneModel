clear all;


ModelParams = struct();
rndSeed = 108;

domainBoundary = struct();
domainBoundary.x_max =  50;
domainBoundary.y_max =  50;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 500;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 16;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 3;

% T-cells
ModelParams.P_A = 1; % Antigen gain rate
ModelParams.P_D = 0.0;  % Antigen loss rate
ModelParams.activatedAge  = 50;

ModelParams.dx_ABM = 0.5;
ModelParams.dy_ABM = 0.5;
ModelParams.dt_ABM = 0.1;            % Time step
ModelParams.da_ABM = 2;            % Time step

ModelParams.NumbAgents = 5000;

ModelParams.dx_PDE = 0.5;
ModelParams.dy_PDE = 0.5;
ModelParams.dt_PDE = 0.1;
ModelParams.da_PDE = 2;

(ModelParams.dx_PDE*ModelParams.dy_PDE)/(4*ModelParams.dt_PDE)

ModelParams.plot_traj = true;
recordRate = 25;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;


rndSeed = 1;

[u, C, dCdx, dCdy, Ic, A_Ic, Inds_Ic_st, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams);

x = linspace(0, params.Lx, params.Nx);
y = linspace(0, params.Ly, params.Ny);
a = linspace(0, params.La, params.Na);

tic;
for n=1:params.nt
     u = computePhenotypeModel(u, dCdx, dCdy, A, Inds_Ic_st, A_Ic, params);
    activation_proportion(n) = squeeze(sum(u,[1,2]))'*a'/params.La;
    
    if true && mod(n,100)==0
        clf;
    
        subplot(2,2,1)
        u_slice = squeeze(sum(u,3));
        surf(x, y, u_slice', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        colormap('jet');
        view(2)
        title('Cell density')
        axis([0 params.Lx 0 params.Ly])
    
        av_a = zeros(params.Nx, params.Ny);
        for ix=1:params.Nx
            for iy=1:params.Ny
                av_a(ix,iy) = a*squeeze(u(ix,iy,:));
            end
        end
        av_a = av_a/(params.La*sum(u(:)));

    
        subplot(2,2,2)
        surf(x, y, av_a', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        colormap('jet');
        view(2)
        axis([0 params.Lx 0 params.Ly])
        title('Antigen density')
        
        
        subplot(2,2,3)
        plot(1:n,activation_proportion(1:n), 'LineWidth',2)
        axis([0 params.nt 0 1])
        title('T-cell activation')

        subplot(2,2,4)
        plot(a, squeeze(sum(u,[1,2])), 'LineWidth',2)
        title('Antigen distribution')
    

        
        sgtitle(['time = ', num2str(params.dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end


end
toc;






