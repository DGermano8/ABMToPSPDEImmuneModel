% 2D Spatial Phenotypic Advection Diffusion Equation Solver (Matlab)

% 1. Define Parameters
Dx = 0.1;       % Spatial diffusion coefficient
chix = 0.05;     % Chemotaxis coefficient
chia = 1;     % Phenotypic advection coefficient

% Define the chemoattractant concentration C(x, y)
% C = @(x, y) exp(-((x - 0.5).^2 + (y - 0.5).^2) / (2 * 0.2^2));
C = @(x, y) 1-x;

% Define the phenotypic advection function f(a)
f = @(a) 1-a; % Example: advection towards a = 0.5

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
t_end = 1.0;
dt = 0.001;
num_steps = floor((t_end - t_start) / dt);

% 2. Initialize Solution
u = zeros(nx, ny, na);
% Initial condition (e.g., a Gaussian distribution in the center)
% center_x_idx = round(nx / 2);
% center_y_idx = round(ny / 2);
% center_a_idx = round(na / 2);
% sigma_x = 5;
% sigma_y = 5;
% sigma_a = 2;
% for i = 1:nx
%     for j = 1:ny
%         for k = 1:na
%             u(i, j, k) = exp(-((i - center_x_idx)^2 / (2 * sigma_x^2) + ...
%                                (j - center_y_idx)^2 / (2 * sigma_y^2) + ...
%                                (k - center_a_idx)^2 / (2 * sigma_a^2)));
%         end
%     end
% end
u(1,1,1) = 1;


% Calculate gradients of C
[gradCx, gradCy] = gradient(C(X, Y), dx, dy);

figure;
tic;
% 3. Time Loop
for n = 1:num_steps
    u_new = u; % Initialize u_new with the current u


    % Loop through each cell
    for i = 1:nx
        for j = 1:ny
            for k = 1:na
                dudt = 0.0;

                % --- Spatial Diffusion Term ---
                laplacian_u = 0.0;
                % x-direction
                flux_west_diff = 0;
                if i > 1
                    flux_west_diff = Dx * (u(i, j, k) - u(i-1, j, k)) / dx;
                elseif i == 1 % No flux boundary
                    flux_west_diff = 0;
                end

                flux_east_diff = 0;
                if i < nx
                    flux_east_diff = Dx * (u(i+1, j, k) - u(i, j, k)) / dx;
                elseif i == nx % No flux boundary
                    flux_east_diff = 0;
                end
                dudt = dudt + (flux_east_diff - flux_west_diff) / dx;

                % y-direction
                flux_south_diff = 0;
                if j > 1
                    flux_south_diff = Dx * (u(i, j, k) - u(i, j-1, k)) / dy;
                elseif j == 1 % No flux boundary
                    flux_south_diff = 0;
                end

                flux_north_diff = 0;
                if j < ny
                    flux_north_diff = Dx * (u(i, j+1, k) - u(i, j, k)) / dy;
                elseif j == ny % No flux boundary
                    flux_north_diff = 0;
                end
                dudt = dudt + (flux_north_diff - flux_south_diff) / dy;

                % --- Chemotaxis Term ---
                div_chemotaxis = 0.0;
                % x-direction
                flux_west_chem = 0;
                if i > 1
                    u_west = 0.5 * (u(i, j, k) + u(i-1, j, k)); % Central difference for u at face
                    flux_west_chem = -chix * u_west * gradCx(j, i-1); % Note the index for gradCx
                elseif i == 1 % No flux boundary
                    flux_west_chem = 0;
                end

                flux_east_chem = 0;
                if i < nx
                    u_east = 0.5 * (u(i+1, j, k) + u(i, j, k)); % Central difference for u at face
                    flux_east_chem = -chix * u_east * gradCx(j, i);   % Note the index for gradCx
                elseif i == nx % No flux boundary
                    flux_east_chem = 0;
                end
                div_chemotaxis = div_chemotaxis + (flux_east_chem - flux_west_chem) / dx;

                % y-direction
                flux_south_chem = 0;
                if j > 1
                    u_south = 0.5 * (u(i, j, k) + u(i, j-1, k)); % Central difference for u at face
                    flux_south_chem = -chix * u_south * gradCy(j-1, i); % Note the index for gradCy
                elseif j == 1 % No flux boundary
                    flux_south_chem = 0;
                end

                flux_north_chem = 0;
                if j < ny
                    u_north = 0.5 * (u(i, j+1, k) + u(i, j, k)); % Central difference for u at face
                    flux_north_chem = -chix * u_north * gradCy(j, i);   % Note the index for gradCy
                elseif j == ny % No flux boundary
                    flux_north_chem = 0;
                end
                div_chemotaxis = div_chemotaxis + (flux_north_chem - flux_south_chem) / dy;

                dudt = dudt - div_chemotaxis;

                % --- Phenotypic Advection Term ---
                adv_phenotype = 0.0;
                f_k_plus_half = f((a(min(k + 1, na)) + a(k)) / 2);
                f_k_minus_half = f((a(k) + a(max(1, k - 1))) / 2);

                u_upwind_plus = 0;
                if chia * f_k_plus_half > 0
                    u_upwind_plus = u(i, j, k);
                else
                    u_upwind_plus = u(i, j, min(k + 1, na));
                end

                u_upwind_minus = 0;
                if chia * f_k_minus_half > 0
                    u_upwind_minus = u(i, j, max(1, k - 1));
                else
                    u_upwind_minus = u(i, j, k);
                end

                flux_a_plus = chia * f_k_plus_half * u_upwind_plus;
                flux_a_minus = chia * f_k_minus_half * u_upwind_minus;

                if k == 1 % Boundary in a (no flux)
                    flux_a_minus = 0;
                end
                if k == na % Boundary in a (no flux)
                    flux_a_plus = 0;
                end

                adv_phenotype = -(flux_a_plus - flux_a_minus) / da;
                dudt = dudt + adv_phenotype;

                % Update u
                u_new(i, j, k) = u(i, j, k) + dt * dudt;
            end
        end
    end
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