% 2D Spatial Phenotypic Advection Diffusion Equation Solver (Matlab)

% 1. Define Parameters
Dx = 0.1;       % Spatial diffusion coefficient
chix = 0.05;     % Chemotaxis coefficient
chia = 0.1;     % Phenotypic advection coefficient


% Define the spatial domain
xmin = 0.0; xmax = 1.0;
ymin = 0.0; ymax = 1.0;
nx = 50;
ny = 50;
dx = (xmax - xmin) / (nx - 1);
dy = (ymax - ymin) / (ny - 1);
x = linspace(xmin, xmax, nx);
y = linspace(ymin, ymax, ny);
[X, Y] = meshgrid(x, y);

% Define the phenotypic domain
amin = 0.0; amax = 1.0;
na = 50;
da = (amax - amin) / (na - 1);
a = linspace(amin, amax, na);

% Define time parameters
t_start = 0.0;
t_end = 10.0;
dt = 0.001;
num_steps = floor((t_end - t_start) / dt);

% 2. Initialize Solution
u = zeros(nx, ny, na);
u(1,:,1) = 1;
u = u/sum(u(:));

Ic = zeros(nx,ny,na);
Ic(round(0.5*nx):end,:,:) = 1;


rho_minus = 0;
rho_plus = 1;

% Define the chemoattractant concentration C(x, y)
[Y,X, A] = meshgrid(y,x,a);
C_m = X/xmax;

figure;
tic;
% 3. Time Loop
for n = 1:num_steps
    u_new = u; % Initialize u_new with the current u


    % Loop through each cell
    % for i = 1:nx
        % for j = 1:ny
            for k = 1:na
                dudt = 0.0;

                % --- Spatial Diffusion Term ---
                % x-direction

                flux_west_diff = Dx * (u(2:end, :, k) - u(1:end-1, :, k)) / dx;
                flux_west_diff = cat(1, zeros(1,ny), flux_west_diff);

                flux_east_diff = Dx * (u(2:end, :, k) - u(1:end-1, :, k)) / dx;
                flux_east_diff = cat(1, flux_east_diff, zeros(1,ny));
                laplacian_u = (flux_east_diff - flux_west_diff) / dx;

                % y-direction
                flux_south_diff = Dx * (u(:, 2:end, k) - u(:, 1:end-1, k)) / dy;
                flux_south_diff = cat(2, zeros(nx,1), flux_south_diff);
                
                flux_north_diff = Dx * (u(:, 2:end, k) - u(:, 1:end-1, k)) / dy;
                flux_north_diff = cat(2, flux_north_diff, zeros(nx,1));
                
                laplacian_u = laplacian_u + (flux_north_diff - flux_south_diff) / dy;

                dudt = dudt + laplacian_u;

                % --- Chemotaxis Term ---
                % x-direction
                dCdx = (C_m(2:end, :, k) - C_m(1:(end-1),:, k)) / dx;

                u_west = 0.5 * (u(2:end, :, k) + u(1:end-1, :, k)); % Central difference for u at face
                flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx
                flux_west_chem = cat(1, zeros(1,ny), flux_west_chem);

                u_east = 0.5 * (u(2:end, :, k) + u(1:end-1, :, k)); % Central difference for u at face
                flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx
                flux_east_chem = cat(1, flux_east_chem, zeros(1,ny));
                
                div_chemotaxis = (flux_east_chem - flux_west_chem) / dx;

                % y-direction
                dCdy = (C_m(:, 2:end,k) - C_m(:,1:(end-1),k)) / dy;

                u_south = 0.5 * (u(:, 2:end, k) + u(:, 1:end-1, k)); % Central difference for u at face
                flux_south_chem = -chix .* u_south .* dCdy; % Note the index for gradCy
                flux_south_chem = cat(2, zeros(nx,1), flux_south_chem);

                u_north = 0.5 * (u(:, 2:end, k) + u(:, 1:end-1, k)); % Central difference for u at face
                flux_north_chem = -chix .* u_north .* dCdy;   % Note the index for gradCy
                flux_north_chem  = cat(2, flux_north_chem, zeros(nx,1));

                div_chemotaxis = div_chemotaxis + (flux_north_chem - flux_south_chem) / dy;

                dudt = dudt - div_chemotaxis;

                % --- Phenotypic Advection Term ---
                a_plus_half_a = (A(:,:,min(k + 1, na)) + A(:,:,k)) / 2;
                f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* Ic(:,:, k) - rho_minus .* a_plus_half_a./amax .* (1-Ic(:,:, k)) ) ;
                % f_k_plus_half = (1 - a_plus_half_a./amax) .*  Ic(:,:, k);
                % f_k_plus_half = f((a(min(k + 1, na)) + a(k)) / 2);
                u_upwind_plus = u(:, :, min(k + 1, na)) .*(chia * f_k_plus_half <= 0) + u(:, :, k) .* (chia * f_k_plus_half > 0);
                
                a_minus_half_a = (A(:,:,k) + A(:,:,(max(1, k - 1))) ) / 2;
                % a_minus_half_a = (a(k) + a((max(1, k - 1))) ) / 2;
                f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* Ic(:,:, k) - rho_minus .* a_minus_half_a./amax .* (1-Ic(:,:, k)) ) ;
                % f_k_minus_half = (1 - a_minus_half_a./amax) .*  Ic(:,:, k);
                % f_k_minus_half = f((a(k) + a(max(1, k - 1))) / 2);
                u_upwind_minus = u(:, :, k) .*(chia * f_k_minus_half <= 0) + u(:, :, max(1, k - 1)) .* (chia * f_k_minus_half > 0);

                flux_a_plus = chia .* f_k_plus_half .* u_upwind_plus;
                flux_a_minus = chia .* f_k_minus_half .* u_upwind_minus;

                if k == 1 % Boundary in a (no flux)
                    flux_a_minus = zeros(nx,ny,1);
                end
                if k == na % Boundary in a (no flux)
                    flux_a_plus = zeros(nx,ny,1);
                end

                adv_phenotype = -(flux_a_plus - flux_a_minus) / da;
                dudt = dudt + adv_phenotype;

                % Update u
                u_new(:, :, k) = u(:, :, k) + dt * dudt;
            end
        % end
    % end
    u = u_new;

    % Optional: Visualize the solution at certain time steps
    if mod(n, 100) == 0
        clf;
    
        subplot(2,2,1)
        u_slice = squeeze(sum(u,3));
        surf(x, y, u_slice', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        colormap('jet');
        view(2)
        axis([0 ymax 0 xmax])
    
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
        colormap('jet');
        view(2)
        axis([0 ymax 0 xmax])
        

        subplot(2,2,4)
        plot(a, squeeze(sum(u,[1,2])))
        axis([0 max(a) 0 max(squeeze(sum(u,[1,2])))])

    
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