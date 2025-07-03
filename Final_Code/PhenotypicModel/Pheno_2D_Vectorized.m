% 2D Spatial Phenotypic Advection Diffusion Equation Solver (Matlab)

% 1. Define Parameters
Dx = 0.5;       % Spatial diffusion coefficient
chix = 1.0;     % Chemotaxis coefficient
chia = 1.0;     % Phenotypic advection coefficient


% Define the spatial domain
xmin = 0.0; xmax = 10.0;
ymin = 0.0; ymax = 1.0;
nx = 51;
ny = 6;
dx = (xmax - xmin) / (nx - 1);
dy = (ymax - ymin) / (ny - 1);
x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
[X, Y] = meshgrid(x, y);

% Define the phenotypic domain
amin = 0.0; amax = 50.0;
na = 101;
da = (amax - amin) / (na - 1);
a = linspace(amin, amax, na);

% Define time parameters
t_start = 0.0;
t_end = 100.0;
dt = 0.005;
num_steps = floor((t_end - t_start) / dt);

% 2. Initialize Solution
u = zeros(nx, ny, na);
u(1,:,1) = 1;
u = u/sum(u(:));

Ic = zeros(nx,ny,na);
Ic(round(0.5*nx):end,:,:) = 1;


rho_minus = 0.9;
rho_plus = 0.5;

% Define the chemoattractant concentration C(x, y)
[Y,X, A] = meshgrid(y,x,a);
C_m = X/xmax;

figure;
tic;
% 3. Time Loop
for n = 1:num_steps
    u_new = u; % Initialize u_new with the current u

        % --- Spatial Diffusion Term ---
        % x-direction

        flux_west_diff = Dx * (u(2:end, :, :) - u(1:end-1, :, :)) / dx;
        flux_west_diff = cat(1, zeros(1,ny,na), flux_west_diff);

        flux_east_diff = Dx * (u(2:end, :, :) - u(1:end-1, :, :)) / dx;
        flux_east_diff = cat(1, flux_east_diff, zeros(1,ny,na));
        laplacian_u = (flux_east_diff - flux_west_diff) / dx;

        % y-direction
        flux_south_diff = Dx * (u(:, 2:end, :) - u(:, 1:end-1, :)) / dy;
        flux_south_diff = cat(2, zeros(nx,1,na), flux_south_diff);
        
        flux_north_diff = Dx * (u(:, 2:end, :) - u(:, 1:end-1, :)) / dy;
        flux_north_diff = cat(2, flux_north_diff, zeros(nx,1,na));
        
        laplacian_u = laplacian_u + (flux_north_diff - flux_south_diff) / dy;


        % --- Chemotaxis Term ---
        % x-direction
        dCdx = (C_m(2:end, :, :) - C_m(1:(end-1),:, :)) / dx;

        u_west = 0.5 * (u(2:end, :, :) + u(1:end-1, :, :)); % Central difference for u at face
        flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx
        flux_west_chem = cat(1, zeros(1,ny,na), flux_west_chem);

        u_east = 0.5 * (u(2:end, :, :) + u(1:end-1, :, :)); % Central difference for u at face
        flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx
        flux_east_chem = cat(1, flux_east_chem, zeros(1,ny,na));
        
        div_chemotaxis = (flux_east_chem - flux_west_chem) / dx;

        % y-direction
        dCdy = (C_m(:, 2:end,:) - C_m(:,1:(end-1),:)) / dy;

        u_south = 0.5 * (u(:, 2:end, :) + u(:, 1:end-1, :)); % Central difference for u at face
        flux_south_chem = -chix .* u_south .* dCdy; % Note the index for gradCy
        flux_south_chem = cat(2, zeros(nx,1,na), flux_south_chem);

        u_north = 0.5 * (u(:, 2:end, :) + u(:, 1:end-1, :)); % Central difference for u at face
        flux_north_chem = -chix .* u_north .* dCdy;   % Note the index for gradCy
        flux_north_chem  = cat(2, flux_north_chem, zeros(nx,1,na));

        div_chemotaxis = div_chemotaxis + (flux_north_chem - flux_south_chem) / dy;


        % --- Phenotypic Advection Term ---
        a_ind_1 = [2:na, na];
        a_plus_half_a = (A(:,:,a_ind_1) + A(:,:,:)) / 2;
        f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* Ic(:,:, :) - rho_minus .* a_plus_half_a./amax .* (1-Ic(:,:, :)) ) ;
        u_upwind_plus = u(:, :, a_ind_1) .*(chia .* f_k_plus_half <= 0) + u(:, :, :) .* (chia .* f_k_plus_half > 0);
        
        a_ind_2 = [1 1:(na-1)];
        a_minus_half_a = (A(:,:,:) + A(:,:,a_ind_2) ) / 2;
        f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* Ic(:,:, :) - rho_minus .* a_minus_half_a./amax .* (1-Ic(:,:, :)) ) ;
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
        view(2)
        title('Cell Density')
        axis([0 xmax 0 ymax])
    
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
        axis([0 xmax 0 ymax])

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