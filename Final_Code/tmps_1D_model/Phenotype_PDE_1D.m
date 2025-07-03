
% Define parameters
Np = 500;

Lx = 10; 
La = 50; 
T  = 25;

Dx = 0.25; 
Da = 0.05; 
Dt = 0.001;

Nt = T/Dt;
Na = La/Da;
Nx = Lx/Dx;

D = 1;
chi = 0.5; 
Da_diffusion = 0.001;
chia = 1; % 'a'-chemotaxis sensitivity
rho_plus = 1;
rho_minus = 0.2;

% Create grid
x = linspace(0, Lx, Nx); 
a = linspace(0, La, Na); 
t = linspace(0, T, Nt + 1);

% Initialize u (example initial condition)
u = zeros(Nx, Na);
u(1,1) = 1;

activation_proportion = zeros(Nt + 1,1);

% Initialize C (example function of x and t)
C = @(xx, tt) xx/Lx;

% Define f(u, a) (example function)
Ic = zeros(Nx,Na);
Ic(round(0.5*Nx):end,:) = 1;
% Ic(:,:) = 1;
f = @(uu, aa, ix, ia)  (rho_plus * (1 - aa/La) * Ic(ix,ia) - rho_minus * aa/La * (1-Ic(ix,ia)) );

% Store solution over time
% u_history = zeros(Nt + 1, Nx, Na);
% u_history(1, :, :) = u;

figure;
% Time stepping loop
for n = 1:Nt
    
    u_new = u;
    t_current = Dt*n;

    for ix = 1:Nx
        for ia = 1:Na
            % Fluxes in x-direction
            flux_x_plus = 0;
            flux_x_minus = 0;
            if ix < Nx
                u_avg_plus_x = (u(ix + 1, ia) + u(ix, ia)) / 2;
                dudx_plus = (u(ix + 1, ia) - u(ix, ia)) / Dx;
                dCdx_plus = (C(x(ix + 1), t_current) - C(x(ix), t_current)) / Dx;
                flux_x_plus = D * dudx_plus - u_avg_plus_x * chi * dCdx_plus;
            end
            if ix > 1
                u_avg_minus_x = (u(ix, ia) + u(ix - 1, ia)) / 2;
                dudx_minus = (u(ix, ia) - u(ix - 1, ia)) / Dx;
                dCdx_minus = (C(x(ix), t_current) - C(x(ix - 1), t_current)) / Dx;
                flux_x_minus = D * dudx_minus - u_avg_minus_x * chi * dCdx_minus;
            end

            % Boundary conditions for x
            if ix == 1
                % flux_x_minus = D * (u(ix + 1, ia) - u(ix, ia)) / Dx - ((u(ix + 1, ia) + u(ix, ia)) / 2) * chi * (C(x(ix + 1), t_current) - C(x(ix), t_current)) / Dx; % Reflective
                flux_x_minus = 0;
            elseif ix == Nx
                % flux_x_plus = D * (u(ix, ia) - u(ix - 1, ia)) / Dx - ((u(ix, ia) + u(ix - 1, ia)) / 2) * chi * (C(x(ix), t_current) - C(x(ix - 1), t_current)) / Dx; % Reflective
                flux_x_plus = 0;
            end

            div_x_final = (flux_x_plus - flux_x_minus) / Dx;



            % Fluxes in a-direction
            flux_a_plus = 0;
            flux_a_minus = 0;
            if ia < Na
                u_avg_plus_a = (u(ix, ia + 1) + u(ix, ia)) / 2;
                a_plus_half_a = a(ia) + Da / 2;
                duda_plus = (u(ix, ia + 1) - u(ix, ia)) / Da;
                flux_a_plus =  chia * f(u_avg_plus_a, a_plus_half_a, ix, ia+1 ) * u_avg_plus_a - Da_diffusion * duda_plus;
            end
            if ia > 1
                u_avg_minus_a = (u(ix, ia) + u(ix, ia - 1)) / 2;
                a_minus_half_a = a(ia) - Da / 2;
                duda_minus = (u(ix, ia) - u(ix, ia - 1)) / Da;
                flux_a_minus = chia * f(u_avg_minus_a, a_minus_half_a, ix, ia-1 ) * u_avg_minus_a - Da_diffusion * duda_minus;
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
            u_new(ix, ia) = u(ix, ia) + Dt * (div_x_final - div_a_final);
        end
    end
    u = u_new;
    % u_history(n + 1, :, :) = u;

    activation_proportion(n) = sum((u*a') / sum(u(:)))/La;

    if mod(n,Np) == 0
        clf;
        
        subplot(3,1,1)
        u_slice = squeeze(u(:, :));
        surf(x, a, u_slice', 'EdgeColor','none');
        % imagesc(x, a, u_slice');

        colorbar;
        xlabel('x');
        ylabel('a');
        colormap('jet');
        view(2)



        subplot(3,2,3)
        plot(a ,sum(u,1), LineWidth=2)
        title('antigen distribution')
        axis([min(a) max(a) 0 1.1*max(sum(u,1))])

        subplot(3,2,4)
        plot(x ,sum(u,2), LineWidth=2)
        title('cell density')
        axis([min(x) max(x) 0 1.1*max(sum(u,2))])

        avg_a = (u*a') / sum(u(:));

        subplot(3,2,5)
        plot(x ,avg_a, LineWidth=2)
        title('average antigen density')
        axis([min(x) max(x) 0 1.1*max(avg_a)])

        subplot(3,2,6)
        plot(1:n, activation_proportion(1:n), LineWidth=2)
        title('activation proportion')
        axis([0 Nt 0 1])

        
        sgtitle(['Time = ', num2str(t_current, '%.2f'), ' mass = ', num2str(sum(u(:)), '%.2f')]);

        drawnow
    end
end
