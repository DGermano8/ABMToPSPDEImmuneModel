% 2D Spatial Phenotypic Advection Diffusion Equation Solver (Matlab)
clear all;
close all;

% 1. Define Parameters
Dx = 0.5;       % Spatial diffusion coefficient
chix = 0.50;     % Chemotaxis coefficient
chia = 1.0;     % Phenotypic advection coefficient

% Define the spatial domain


% Define the phenotypic domain
amin = 0.0; amax = 20.0;
da =  0.2;
na = (amax - amin)/da + 1;
a = linspace(amin, amax, na);

rho_minus = 0.0;
rho_plus = 1;


ModelParams = struct();
rndSeed = 108;

domainBoundary = struct();
domainBoundary.x_max =  7;
domainBoundary.y_max =  7;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 200;
ModelParams.p_move = Dx;

% ChemoTaxis params
ModelParams.C_sens = chix;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 1;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = rho_plus; % Antigen gain rate
ModelParams.P_D = rho_minus;  % Antigen loss rate
ModelParams.activatedAge  = amax;


dx = 0.25;
dy = 0.25;
dt = 0.01;

ModelParams.dx_PDE = dx;
ModelParams.dy_PDE = dx;
ModelParams.dt_PDE = dt;
ModelParams.da_PDE = da;

nx = (domainBoundary.x_max-domainBoundary.x_min)/ModelParams.dx_PDE + 1;
ny = (domainBoundary.y_max-domainBoundary.y_min)/ModelParams.dy_PDE + 1;
x = linspace(domainBoundary.x_min, domainBoundary.x_max, nx);
y = linspace(domainBoundary.y_min, domainBoundary.y_max, ny);

% Define the chemoattractant concentration C(x, y)
[Y,X, A] = meshgrid(y,x,a);

[u, C_m, dCdx, dCdy, Ic, A_Ic, Inds_Ic, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams);

a_ind_1 = params.a_ind_1;
a_ind_2 = params.a_ind_2;

inds = Inds_Ic.inds;
inds_north = Inds_Ic.inds_north;
inds_south = Inds_Ic.inds_south;
inds_east = Inds_Ic.inds_east;
inds_west = Inds_Ic.inds_west;


%%
figure;
tic;
% 3. Time Loop
for n = 1:params.nt

    % --- Spatial Flux Term ---
    % x-direction
    flux_west_diff = Dx * (u(2:end, :, :) - u(1:end-1, :, :)) / dx;
    u_west = 0.5 * (u(2:end, :, :) + u(1:end-1, :, :)); % Central difference for u at face
    flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx

    flux_west = flux_west_chem + flux_west_diff;
    flux_west = cat(1, zeros(1,ny,na), flux_west);
    flux_west(inds_west) = 0;


    flux_east_diff = Dx * (u(2:end, :, :) - u(1:end-1, :, :)) / dx;
    u_east = 0.5 * (u(2:end, :, :) + u(1:end-1, :, :)); % Central difference for u at face
    flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx
    
    flux_east = flux_east_chem + flux_east_diff;
    flux_east = cat(1, flux_east, zeros(1,ny,na));
    flux_east(inds_east) = 0;

    grad_x = (flux_east - flux_west) / dx;


    % y-direction
    flux_south_diff = Dx * (u(:, 2:end, :) - u(:, 1:end-1, :)) / dy;
    u_south = 0.5 * (u(:, 2:end, :) + u(:, 1:end-1, :)); % Central difference for u at face
    flux_south_chem = -chix .* u_south .* dCdy; % Note the index for gradCy
    
    flux_south = flux_south_chem + flux_south_diff;
    flux_south = cat(2, zeros(nx,1,na), flux_south);
    flux_south(inds_south) = 0;


    flux_north_diff = Dx * (u(:, 2:end, :) - u(:, 1:end-1, :)) / dy;
    u_north = 0.5 * (u(:, 2:end, :) + u(:, 1:end-1, :)); % Central difference for u at face
    flux_north_chem = -chix .* u_north .* dCdy;   % Note the index for gradCy

    flux_north = flux_north_chem + flux_north_diff;
    flux_north  = cat(2, flux_north, zeros(nx,1,na));
    flux_north(inds_north) = 0;

    grad_y = (flux_north - flux_south) / dy;


    % --- Phenotypic Advection Term ---
    a_plus_half_a = (A(:,:,a_ind_1) + A(:,:,:)) / 2;
    f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_plus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
    u_upwind_plus = u(:, :, a_ind_1) .*(f_k_plus_half <= 0) + u(:, :, :) .* (f_k_plus_half > 0);
    
    a_minus_half_a = (A(:,:,:) + A(:,:,a_ind_2) ) / 2;
    f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_minus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
    u_upwind_minus = u(:, :, :) .*(f_k_minus_half <= 0) + u(:, :, a_ind_2) .* (f_k_minus_half > 0);

    flux_a_plus = chia .* f_k_plus_half .* u_upwind_plus;
    flux_a_minus = chia .* f_k_minus_half .* u_upwind_minus;

    % Boundary in a (no flux)
    flux_a_minus(:,:,1) = 0;
    % Boundary in a (no flux)
    flux_a_plus(:,:,end) = 0;

    grad_a = -(flux_a_plus - flux_a_minus) / da;


    % Update u
    u = u + dt * (grad_x + grad_y + grad_a);
    u(inds) = 0;

    % Optional: Visualize the solution at certain time steps
    if mod(n, 500) == 0
        clf;
    
        subplot(2,2,1)
        u_slice = squeeze(sum(u,3));
        surf(x, y, u_slice', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        % colormap('jet');
        % hold on; scatter3(Inds_Ic(:,1)+1,Inds_Ic(:,2)+1,ones(length(Inds_Ic(:,1))),'magenta')
        view(2)
        title('Cell Density')
        axis([0 domainBoundary.x_max 0 domainBoundary.y_max])
    
        av_a = zeros(nx, ny);
        for ix=1:nx
            for iy=1:ny
                av_a(ix,iy) = a*squeeze(u(ix,iy,:));
            end
        end
        av_a = av_a/amax*sum(u(:));

    
        subplot(2,2,2)
        surf(x, y, av_a', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        % colormap('jet');
        view(2)
        title('Average Antigen Level')
        axis([0 domainBoundary.x_max 0 domainBoundary.y_max])

        subplot(2,2,3)
        surf(x, y, C_m(:,:,1)', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        view(2)
        
        subplot(2,2,4)
        plot(a, squeeze(sum(u,[1,2])), 'LineWidth',2)
        axis([0 max(a) 0 max(squeeze(sum(u,[1,2])))])
        title('Distribution of Cell Antigen')
    
        sgtitle(['time = ', num2str(dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end
end
toc;

disp('Simulation finished.');

%%
clf;

subplot(2,2,1)
u_slice = squeeze(sum(u,3));
surf(x, y, u_slice', 'EdgeColor','none');
colorbar;
xlabel('x');
ylabel('y');
% colormap('jet');
% hold on; scatter3(Inds_Ic(:,1)+1,Inds_Ic(:,2)+1,ones(length(Inds_Ic(:,1))),'magenta')
view(2)
title('Cell Density')
axis([0 domainBoundary.x_max 0 domainBoundary.y_max])

av_a = zeros(nx, ny);
for ix=1:nx
    for iy=1:ny
        av_a(ix,iy) = a*squeeze(u(ix,iy,:));
    end
end
av_a = av_a/amax*sum(u(:));


subplot(2,2,2)
surf(x, y, av_a', 'EdgeColor','none');
colorbar;
xlabel('x');
ylabel('y');
% colormap('jet');
view(2)
title('Average Antigen Level')
axis([0 domainBoundary.x_max 0 domainBoundary.y_max])

subplot(2,2,3)
surf(x, y, C_m(:,:,1)', 'EdgeColor','none');
colorbar;
xlabel('x');
ylabel('y');
view(2)

subplot(2,2,4)
plot(a, squeeze(sum(u,[1,2])), 'LineWidth',2)
axis([0 max(a) 0 max(squeeze(sum(u,[1,2])))])
title('Distribution of Cell Antigen')

sgtitle(['time = ', num2str(dt*n), ' mass = ', num2str(sum(u(:)))])
drawnow;