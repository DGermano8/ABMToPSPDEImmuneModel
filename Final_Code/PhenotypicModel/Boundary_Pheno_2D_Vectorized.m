% 2D Spatial Phenotypic Advection Diffusion Equation Solver (Matlab)
clear all;
close all;

% 1. Define Parameters
Dx = 0.5;       % Spatial diffusion coefficient
chix = 0.50;     % Chemotaxis coefficient
chia = 1.0;     % Phenotypic advection coefficient

% Define the spatial domain


% Define the phenotypic domain
amin = 0.0; amax = 50.0;
na = 101;
da = (amax - amin) / (na - 1);
a = linspace(amin, amax, na);

rho_minus = 0.5;
rho_plus = 1;




ModelParams = struct();
rndSeed = 108;

domainBoundary = struct();
domainBoundary.x_max =  50;
domainBoundary.y_max =  50;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 100;
ModelParams.p_move = Dx;

% ChemoTaxis params
ModelParams.C_sens = chix;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 8;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = rho_plus; % Antigen gain rate
ModelParams.P_D = rho_minus;  % Antigen loss rate
ModelParams.activatedAge  = amax;



dx = 0.2;
dy = 0.2;
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


[u, C_m, Ic, A_Ic, Inds_Ic, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams);

% Inds_Ic = [];
% Inds_Ic = [10*ones(nx,1)'; 1:nx']';

inds = [];
inds_north = [];
inds_south = [];
inds_east = [];
inds_west = [];
for ii = 1:size(Inds_Ic)
    for ia = 1:na
        
        if (Ic(Inds_Ic(ii,1)+1,Inds_Ic(ii,2)) + Ic(Inds_Ic(ii,1)-1,Inds_Ic(ii,2)) ...
                + Ic(Inds_Ic(ii,1),Inds_Ic(ii,2)+1) + Ic(Inds_Ic(ii,1),Inds_Ic(ii,2)-1)) > 0
            inds = [inds; sub2ind(size(u), Inds_Ic(ii,1),Inds_Ic(ii,2),ia)];
        end

        if ~(ismember([Inds_Ic(ii,1),Inds_Ic(ii,2)-1],Inds_Ic,'rows'))
            inds_north = [inds_north; sub2ind(size(u), Inds_Ic(ii,1),Inds_Ic(ii,2)-1,ia)];
        end
        if ~(ismember([Inds_Ic(ii,1),Inds_Ic(ii,2)+1],Inds_Ic,'rows'))
            inds_south = [inds_south; sub2ind(size(u), Inds_Ic(ii,1),Inds_Ic(ii,2)+1,ia)];
        end
        if ~(ismember([Inds_Ic(ii,1)-1,Inds_Ic(ii,2)],Inds_Ic,'rows'))
            inds_east = [inds_east; sub2ind(size(u), Inds_Ic(ii,1)-1,Inds_Ic(ii,2),ia)];
        end
        if ~(ismember([Inds_Ic(ii,1)+1,Inds_Ic(ii,2)],Inds_Ic,'rows'))
            inds_west = [inds_west; sub2ind(size(u), Inds_Ic(ii,1)+1,Inds_Ic(ii,2),ia)];
        end       
    end
end


%%
figure;
tic;
% 3. Time Loop
for n = 1:params.nt
    n
    u_new = u; % Initialize u_new with the current u

        % --- Spatial Diffusion Term ---
        % x-direction

        flux_west_diff = Dx * (u(2:end, :, :) - u(1:end-1, :, :)) / dx;
        flux_west_diff = cat(1, zeros(1,ny,na), flux_west_diff);
        flux_west_diff(inds_west) = 0;

        flux_east_diff = Dx * (u(2:end, :, :) - u(1:end-1, :, :)) / dx;
        flux_east_diff = cat(1, flux_east_diff, zeros(1,ny,na));
        flux_east_diff(inds_east) = 0;

        laplacian_u = (flux_east_diff - flux_west_diff) / dx;

        % y-direction
        flux_south_diff = Dx * (u(:, 2:end, :) - u(:, 1:end-1, :)) / dy;
        flux_south_diff = cat(2, zeros(nx,1,na), flux_south_diff);
        flux_south_diff(inds_south) = 0;

        flux_north_diff = Dx * (u(:, 2:end, :) - u(:, 1:end-1, :)) / dy;
        flux_north_diff = cat(2, flux_north_diff, zeros(nx,1,na));
        flux_north_diff(inds_north) = 0;

        laplacian_u = laplacian_u + (flux_north_diff - flux_south_diff) / dy;


        % --- Chemotaxis Term ---
        % x-direction
        dCdx = (C_m(2:end, :, :) - C_m(1:(end-1),:, :)) / dx;

        u_west = 0.5 * (u(2:end, :, :) + u(1:end-1, :, :)); % Central difference for u at face
        flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx
        flux_west_chem = cat(1, zeros(1,ny,na), flux_west_chem);
        flux_west_chem(inds_west) = 0;

        u_east = 0.5 * (u(2:end, :, :) + u(1:end-1, :, :)); % Central difference for u at face
        flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx
        flux_east_chem = cat(1, flux_east_chem, zeros(1,ny,na));
        flux_east_chem(inds_east) = 0;

        div_chemotaxis = (flux_east_chem - flux_west_chem) / dx;

        % y-direction
        dCdy = (C_m(:, 2:end,:) - C_m(:,1:(end-1),:)) / dy;

        u_south = 0.5 * (u(:, 2:end, :) + u(:, 1:end-1, :)); % Central difference for u at face
        flux_south_chem = -chix .* u_south .* dCdy; % Note the index for gradCy
        flux_south_chem = cat(2, zeros(nx,1,na), flux_south_chem);
        flux_south_chem(inds_south) = 0;

        u_north = 0.5 * (u(:, 2:end, :) + u(:, 1:end-1, :)); % Central difference for u at face
        flux_north_chem = -chix .* u_north .* dCdy;   % Note the index for gradCy
        flux_north_chem  = cat(2, flux_north_chem, zeros(nx,1,na));
        flux_north_chem(inds_north) = 0;

        div_chemotaxis = div_chemotaxis + (flux_north_chem - flux_south_chem) / dy;


        % --- Phenotypic Advection Term ---
        a_ind_1 = [2:na, na];
        a_plus_half_a = (A(:,:,a_ind_1) + A(:,:,:)) / 2;
        f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_plus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
        u_upwind_plus = u(:, :, a_ind_1) .*(chia .* f_k_plus_half <= 0) + u(:, :, :) .* (chia .* f_k_plus_half > 0);
        
        a_ind_2 = [1 1:(na-1)];
        a_minus_half_a = (A(:,:,:) + A(:,:,a_ind_2) ) / 2;
        f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_minus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
        u_upwind_minus = u(:, :, :) .*(chia .* f_k_minus_half <= 0) + u(:, :, a_ind_2) .* (chia .* f_k_minus_half > 0);

        flux_a_plus = chia .* f_k_plus_half .* u_upwind_plus;
        flux_a_minus = chia .* f_k_minus_half .* u_upwind_minus;

        % Boundary in a (no flux)
        flux_a_minus(:,:,1) = 0;
        % Boundary in a (no flux)
        flux_a_plus(:,:,end) = 0;

        adv_phenotype = -(flux_a_plus - flux_a_minus) / da;

        
    % Update u
    u = u + dt * (laplacian_u + div_chemotaxis + adv_phenotype);
    % for ii = 1:length(Inds_Ic)
    %     u(Inds_Ic(ii,1),Inds_Ic(ii,2),:) = 0;
    % end
    u(inds) = 0;

    % Optional: Visualize the solution at certain time steps
    if mod(n, 500) == 0
        clf;
    
        subplot(2,2,1)
        u_slice = squeeze(sum(u,3));
        surf(1:nx, 1:ny, u_slice', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        % colormap('jet');
        % hold on; scatter3(Inds_Ic(:,1)+1,Inds_Ic(:,2)+1,ones(length(Inds_Ic(:,1))),'magenta')
        view(2)
        title('Cell Density')
        % axis([0 domainBoundary.x_max 0 domainBoundary.y_max])
    
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
        plot(a, squeeze(sum(u,[1,2])))
        axis([0 max(a) 0 max(squeeze(sum(u,[1,2])))])
        title('Distribution of Cell Antigen')
    
        sgtitle(['time = ', num2str(dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end
end
toc;

disp('Simulation finished.');

% Optional: Visualize the final solution
% figure(2);
% slice(X, Y, a, u(:, :, :), [0.5], [0.5], [amin, (amin+amax)/2, amax]);
% xlabel('x'); ylabel('y'); zlabel('a');
% title(['Final Solution at Time: ', num2str(t_end)]);
% colorbar;