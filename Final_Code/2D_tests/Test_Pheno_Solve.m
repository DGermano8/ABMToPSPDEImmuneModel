clc;
close all;
clear all;

profile on;

ModelParams = struct();

domainBoundary = struct();
domainBoundary.x_max =  20;
domainBoundary.y_max =  20;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 200;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 2;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 1.0; % Antigen gain rate
ModelParams.P_D = 0.5;  % Antigen loss rate
ModelParams.activatedAge  = 50;

dx = 0.5;
dy = dx;
dt = 0.1;
da = 1;


ModelParams.dx_PDE = dx;
ModelParams.dy_PDE = dy;
ModelParams.dt_PDE = dt;
ModelParams.da_PDE = da;

ModelParams.plot_traj = true;
recordRate = 200;
ModelParams.t_plot = recordRate/ModelParams.dt_PDE;


rndSeed = 2;

[u, ~, dCdx, dCdy, Ic, A_Ic, Inds_Ic, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams);
activation_proportion_pheno = zeros(params.nt,1);
x = linspace(0, params.Lx, params.Nx);
y = linspace(0, params.Ly, params.Ny);
% a = linspace(0, 1, params.Na);
a = linspace(0, params.La, params.Na);


%%

tic;

figure;
for n=0:params.nt
    n/params.nt;
     % u = computePhenotypeModel(u, dCdx, dCdy, A, Inds_Ic, A_Ic, params);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (u, dCdx, dCdy, A, Inds_Ic, A_Ic, params)

    Dx = params.D;
    chix = params.C_chi;
    chia = params.chia;
    rho_plus = params.P_A;
    rho_minus = params.P_D;
    amax = params.activatedAge;

    dx = params.dx;
    dy = params.dy;
    da = params.da;
    dt = params.dt;

    nx = params.Nx;
    ny = params.Ny;
    na = params.Na;

    inds_west = Inds_Ic.inds_west;
    inds_east = Inds_Ic.inds_east;
    inds_south = Inds_Ic.inds_south;
    inds_north = Inds_Ic.inds_north;
    inds = Inds_Ic.inds;

    a_ind_1 = params.a_ind_1;
    a_ind_2 = params.a_ind_2;

    % --- Spatial Flux Term ---
    uxdiff = u(2:end, :, :) - u(1:end-1, :, :);
    uxadd = u(2:end, :, :) + u(1:end-1, :, :);

    % x-direction
    flux_west_diff = Dx * (uxdiff) / dx;
    u_west = 0.5 * (uxadd);
    flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx

    flux_west = flux_west_chem + flux_west_diff;
    flux_west = cat(1, zeros(1,ny,na), flux_west);
    flux_west(inds_west) = 0;


    flux_east_diff = Dx * (uxdiff) / dx;
    u_east = 0.5 * (uxadd); % Central difference for u at face
    flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx

    flux_east = flux_east_chem + flux_east_diff;
    flux_east = cat(1, flux_east, zeros(1,ny,na));
    flux_east(inds_east) = 0;

    grad_x = (flux_east - flux_west) / dx;


    % y-direction
    uydiff = u(:, 2:end, :) - u(:, 1:end-1, :);
    uyadd = u(:, 2:end, :) + u(:, 1:end-1, :);

    flux_south_diff = Dx * (uydiff) / dy;
    u_south = 0.5 * (uyadd); % Central difference for u at face
    flux_south_chem = -chix .* u_south .* dCdy; % Note the index for gradCy

    flux_south = flux_south_chem + flux_south_diff;
    flux_south = cat(2, zeros(nx,1,na), flux_south);
    flux_south(inds_south) = 0;


    flux_north_diff = Dx * (uydiff) / dy;
    u_north = 0.5 * (uyadd); % Central difference for u at face
    flux_north_chem = -chix .* u_north .* dCdy;   % Note the index for gradCy

    flux_north = flux_north_chem + flux_north_diff;
    flux_north  = cat(2, flux_north, zeros(nx,1,na));
    flux_north(inds_north) = 0;

    grad_y = (flux_north - flux_south) / dy;

    % --- Phenotypic Advection Term ---
    a_plus_half_a = (A(:,:,a_ind_1) + A(:,:,:)) / 2;
    f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_plus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
    f_k_plus_half_lessThanZero = (f_k_plus_half <= 0);
    u_upwind_plus = u(:, :, a_ind_1) .*(f_k_plus_half_lessThanZero) + u(:, :, :) .* (1-f_k_plus_half_lessThanZero);

    a_minus_half_a = (A(:,:,:) + A(:,:,a_ind_2) ) / 2;
    f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_minus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
    f_k_minus_half_lessThanZero = (f_k_minus_half <= 0);
    u_upwind_minus = u(:, :, :) .*(f_k_minus_half_lessThanZero) + u(:, :, a_ind_2) .* (1-f_k_minus_half_lessThanZero);

    flux_a_plus = chia .* f_k_plus_half .* u_upwind_plus;
    flux_a_minus = chia .* f_k_minus_half .* u_upwind_minus;

    % Boundary in a (no flux)
    flux_a_minus(:,:,1) = 0;
    % Boundary in a (no flux)
    flux_a_plus(:,:,end) = 0;

    grad_a = -( (flux_a_plus - flux_a_minus) ) / da;

    % Update u
    u = u + dt * (grad_x + grad_y + grad_a);
    u(inds) = 0;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    activation_proportion_pheno(n+1) = squeeze(sum(u,[1,2]))'*a'/params.activatedAge;
       
    if true && mod(n,recordRate)==0
        clf;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ax1 = subplot(1,3,1);
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
    
        ax2 = subplot(1,3,2);
        surf(x, y, av_a', 'EdgeColor','none');
        colorbar;
        view(2)
        ax2.XTick = [];
        ax2.YTick = [];
        axis([0 params.Lx 0 params.Ly])
        xlabel('x', 'FontSize',14, 'Interpreter','latex');
        ylabel('y', 'FontSize',14, 'Interpreter','latex');
        title('PS-PDE antigen density', 'FontSize',14, 'Interpreter','latex');        
        
        subplot(1,3,3)
        plot(ModelParams.dt_PDE*(0:n),activation_proportion_pheno(1:(n+1)), 'LineWidth',2)
        axis([0 params.nt 0 1])
        
        ax3 = subplot(1,3,3);
        hold on
        xlabel('Time', 'FontSize',14, 'Interpreter','latex'); 
        axis([0 ModelParams.T_final 0 1])
        title('Proportion of activation', 'FontSize',14, 'Interpreter','latex')
        
        sgtitle(['time = ', num2str(params.dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end


end

toc;

profile viewer
profile off;

