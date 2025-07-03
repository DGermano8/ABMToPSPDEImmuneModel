% 2D Spatial Phenotypic Advection Diffusion Equation Solver (Matlab)

% 1. Define Parameters
Dx = 0.5;       % Spatial diffusion coefficient
chix = 1.0;     % Chemotaxis coefficient
chia = 1.0;     % Phenotypic advection coefficient


% Define the spatial domain
xmin = 0.0; xmax = 10.0;
dx = 0.5;
nx = xmax/dx + 1;
x = linspace(xmin, xmax, nx);

% Define the phenotypic domain
amin = 0.0; amax = 20.0;

da = 0.1;
na = amax/da + 1;
a = linspace(amin, amax, na);

% Define time parameters
t_start = 0.0;
t_end = 500.0;
dt = 0.02;
num_steps = floor((t_end - t_start) / dt);

% 2. Initialize Solution
u = zeros(nx, na);
u(1,1) = 1;
u = u/sum(u(:));

Ic = zeros(nx,na);
Ic(round(0.5*nx):end,:,:) = 1;


rho_minus = 0.5;
rho_plus = 1;

% Define the chemoattractant concentration C(x, y)
[X, A] = meshgrid(x,a);
X=X'; A=A';
C_m = 1-X/xmax;

figure;
tic;
% 3. Time Loop
for n = 1:num_steps

    % --- Spatial Diffusion Term ---
    % x-direction

    flux_west_diff = Dx * (u(2:end, :) - u(1:end-1, :)) / dx;
    flux_west_diff = cat(1, zeros(1,na), flux_west_diff);

    flux_east_diff = Dx * (u(2:end, :) - u(1:end-1, :)) / dx;
    flux_east_diff = cat(1, flux_east_diff, zeros(1,na));
    laplacian_u = (flux_east_diff - flux_west_diff) / dx;

    % --- Chemotaxis Term ---
    % x-direction
    dCdx = (C_m(2:end, :) - C_m(1:(end-1), :)) / dx;

    u_west = 0.5 * (u(2:end, :) + u(1:end-1, :)); % Central difference for u at face
    flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx
    flux_west_chem = cat(1, zeros(1,na), flux_west_chem);

    u_east = 0.5 * (u(2:end, :) + u(1:end-1, :)); % Central difference for u at face
    flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx
    flux_east_chem = cat(1, flux_east_chem, zeros(1,na));
    
    div_chemotaxis = (flux_east_chem - flux_west_chem) / dx;

    % --- Phenotypic Advection Term ---
    a_ind_1 = [2:na, na];
    a_ind_1 = a_ind_1(:);
    a_plus_half_a = (A(:,a_ind_1) + A(:,:)) / 2;
    f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* Ic(:, :) - rho_minus .* a_plus_half_a./amax .* (1-Ic(:, :)) ) ;
    u_upwind_plus =  u(:, a_ind_1) .*(f_k_plus_half <= 0) + u(:, :) .* (f_k_plus_half > 0) ;
    
    a_ind_2 = [1 1:(na-1)];
    a_ind_2 = a_ind_2(:);
    a_minus_half_a = (A(:,:) + A(:,a_ind_2) ) / 2;
    f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* Ic(:, :) - rho_minus .* a_minus_half_a./amax .* (1-Ic(:, :)) ) ;
    u_upwind_minus = u(:, :) .*( f_k_minus_half <= 0) + u(:, a_ind_2) .* (f_k_minus_half > 0);

    flux_a_plus = chia .* f_k_plus_half .* u_upwind_plus;
    flux_a_minus = chia .* f_k_minus_half .* u_upwind_minus;

    % Boundary in a (no flux)
    flux_a_minus(:,1) = 0;
    % Boundary in a (no flux)
    flux_a_plus(:,end) = 0;

    adv_phenotype = -(flux_a_plus - flux_a_minus) / da;

    u = u + dt * (laplacian_u - div_chemotaxis + adv_phenotype);

    % Optional: Visualize the solution at certain time steps
    if mod(n, 100) == 0
        clf;
    
        subplot(1,3,1)
        u_slice = squeeze(sum(u,2));
        plot(x,  u_slice', 'LineWidth',2);
        xlabel('x');
        title('Cell Density')

        av_a = u*a';  
        subplot(1,3,2)
        plot(x, av_a', 'LineWidth',2);
        xlabel('x');
        title('Average Antigen Level')
        
        subplot(1,3,3)
        plot(a, squeeze(sum(u,1)), 'linewidth',2)
        xlabel('a')
        title('Distribution of Cell Antigen')
    
        sgtitle(['time = ', num2str(dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end
end
toc;

disp('Simulation finished.');
