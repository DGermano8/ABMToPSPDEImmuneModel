
% Define parameters
Np = 500;

Lx = 10; 
Ly = 10; 
La = 10; 
T  = 25;

Dx = 0.25; 
Dy = 0.25; 
Da = 0.1; 
Dt = 0.005;

Nt = T/Dt;
Na = La/Da;
Nx = Lx/Dx;
Ny = Ly/Dy;

D = 0.5;
chi = 0.5; 
Da_diffusion = 0.05;
chia = 1; % 'a'-chemotaxis sensitivity
rho_plus = 1;
rho_minus = 0.2;

% Create grid
x = linspace(0, Lx, Nx); 
y = linspace(0, Ly, Ny); 
a = linspace(0, La, Na); 
t = linspace(0, T, Nt + 1);

% Initialize u (example initial condition)
u = zeros(Nx, Ny, Na);
u(1,1,1) = 1;

activation_proportion = zeros(Nt + 1,1);

% Initialize C (example function of x and t)
C = @(xx, yy, tt) xx/Lx;

% Define f(u, a) (example function)
Ic = zeros(Nx,Ny,Na);
Ic(round(0.5*Nx):end,:,:) = 1;
% Ic(:,:) = 1;
f = @(uu, aa, ix, iy, ia)  (rho_plus * (1 - aa/La) * Ic(ix,iy,ia) - rho_minus * aa/La * (1-Ic(ix,iy,ia)) );

% Store solution over time
% u_history = zeros(Nt + 1, Nx, Na);
% u_history(1, :, :) = u;

figure;
% Time stepping loop
for n = 1:Nt
    
    u_new = u;
    t_current = Dt*n;

    for ix = 1:Nx
        for iy = 1:Ny
            for ia = 1:Na
                % Fluxes in x-direction
                flux_x_plus = 0;
                flux_x_minus = 0;
                if ix < Nx
                    u_avg_plus_x = (u(ix + 1, iy, ia) + u(ix,iy,ia)) / 2;
                    dudx_plus = (u(ix + 1, iy,ia) - u(ix,iy,ia)) / Dx;
                    dCdx_plus = (C(x(ix + 1), y(iy), t_current) - C(x(ix), y(iy), t_current)) / Dx;
                    flux_x_plus = D * dudx_plus - u_avg_plus_x * chi * dCdx_plus;
                end
                if ix > 1
                    u_avg_minus_x = (u(ix,iy, ia) + u(ix - 1,iy, ia)) / 2;
                    dudx_minus = (u(ix,iy, ia) - u(ix - 1,iy, ia)) / Dx;
                    dCdx_minus = (C(x(ix), y(iy), t_current) - C(x(ix - 1), y(iy), t_current)) / Dx;
                    flux_x_minus = D * dudx_minus - u_avg_minus_x * chi * dCdx_minus;
                end
                % Boundary conditions for x
                if ix == 1
                    flux_x_minus = 0;
                elseif ix == Nx
                    flux_x_plus = 0;
                end
                div_x_final = (flux_x_plus - flux_x_minus) / Dx;

                % Fluxes in y-direction
                flux_y_plus = 0;
                flux_y_minus = 0;
                if iy < Ny
                    u_avg_plus_y = (u(ix, iy + 1, ia) + u(ix,iy,ia)) / 2;
                    dudy_plus = (u(ix, iy + 1,ia) - u(ix,iy,ia)) / Dy;
                    dCdy_plus = (C(x(ix), y(iy + 1), t_current) - C(x(ix), y(iy), t_current)) / Dy;
                    flux_y_plus = D * dudy_plus - u_avg_plus_y * chi * dCdy_plus;
                end
                if iy > 1
                    u_avg_minus_y = (u(ix,iy, ia) + u(ix,iy - 1, ia)) / 2;
                    dudy_minus = (u(ix,iy, ia) - u(ix,iy - 1, ia)) / Dy;
                    dCdy_minus = (C(x(ix), y(iy), t_current) - C(x(ix), y(iy - 1), t_current)) / Dy;
                    flux_y_minus = D * dudy_minus - u_avg_minus_y * chi * dCdy_minus;
                end
                % Boundary conditions for x
                if iy == 1
                    flux_y_minus = 0;
                elseif iy == Ny
                    flux_y_plus = 0;
                end
                div_y_final = (flux_y_plus - flux_y_minus) / Dy;
    
    
    
                % Fluxes in a-direction
                flux_a_plus = 0;
                flux_a_minus = 0;
                if ia < Na
                    u_avg_plus_a = (u(ix, iy, ia + 1) + u(ix, iy, ia)) / 2;
                    a_plus_half_a = a(ia) + Da / 2;
                    duda_plus = (u(ix, iy, ia + 1) - u(ix, iy, ia)) / Da;
                    flux_a_plus =  chia * f(u_avg_plus_a, a_plus_half_a, ix, iy, ia+1 ) * u_avg_plus_a - Da_diffusion * duda_plus;
                end
                if ia > 1
                    u_avg_minus_a = (u(ix, iy, ia) + u(ix, iy, ia - 1)) / 2;
                    a_minus_half_a = a(ia) - Da / 2;
                    duda_minus = (u(ix, iy, ia) - u(ix, iy, ia - 1)) / Da;
                    flux_a_minus = chia * f(u_avg_minus_a, a_minus_half_a, ix, iy, ia-1 ) * u_avg_minus_a - Da_diffusion * duda_minus;
                end
    
                % Boundary conditions for a
                if ia == 1
                    % Outflow: flux at a=0 depends on the gradient into the domain
                    % flux_a_minus = chia * f((u(ix, ia + 1) + u(ix, ia)) / 2, a(ia) + Da / 2) * (u(ix, ia + 1) + u(ix, ia)) / 2 - Da_diffusion * (u(ix, ia + 1) - u(ix, ia)) / Da;
                    flux_a_minus = 0;
                elseif ia == Na
                    % Absorbing: flux at a=La uses the value at the boundary
                    % flux_a_plus = chia * f(u(ix, ia), a(ia) + Da / 2) * u(ix, ia) - Da_diffusion * (0 - u(ix, ia)) / Da; % Assume u outside is 0
                    flux_a_plus = 0;
                end
    
                div_a_final = (flux_a_plus - flux_a_minus) / Da;
    
                % Update u
                u_new(ix, iy, ia) = u(ix, iy, ia) + Dt * (div_x_final + div_y_final - div_a_final);
            end
        end
    end
    u = u_new;
    % u_history(n + 1, :, :) = u;

    % activation_proportion(n) = sum((u*a') / sum(u(:)))/La;

    if true && mod(n,10)==0
        clf;
    
        subplot(1,2,1)
        u_slice = squeeze(sum(u,3));
        surf(x, y, u_slice', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        colormap('jet');
        view(2)
    
        av_a = zeros(Nx, Ny);
        for ix=1:Nx
            for iy=1:Ny
                av_a(ix,iy) = a*squeeze(u(ix,iy,:));
            end
        end
        av_a = av_a/La*sum(u(:));
    
        subplot(1,2,2)
        surf(x, y, av_a', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        colormap('jet');
        view(2)
        
    
    
        sgtitle(['time = ', num2str(Dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end
end
