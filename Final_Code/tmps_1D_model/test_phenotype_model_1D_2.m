%% 
% \frac{\partial u}{\partial t} = + D \nabla^2 u  - \nabla \dot (u \chi \nabla C) 
%                               - \chi_a \nabla_a (f(u,a) u)  + D_a \nabla^2_a u

% Define parameters
Nx = 200;          % Number of grid points in x
Na = 100;           % Number of grid points in a
Lx = 100;           % Length of spatial domain
La = 5;            % Length of 'a' domain
T = 10;            % Total time
Dx = Lx / (Nx - 1); % Spatial step size
Da = La / (Na - 1); % 'a' step size
Dt = 0.001;       % Time step size
Nt = T/Dt;
D_x = 1;           % Diffusion coefficient
chi_x = 100;        % Chemotaxis sensitivity

D_a = 0.1;           % Diffusion coefficient
chi_a = 0.1;        % Chemotaxis sensitivity

% f = @(uu, aa) 1 + 0.1 * uu .* aa;
% f = @(uu, aa) 0.1 * uu .* (1 - aa/La);
f = @(uu, aa) 0.1 ;

% Create grid
x = linspace(0, Lx, Nx);
a = linspace(0, La, Na);
t = linspace(0, T, Nt + 1);

% Initialize u (example initial condition)
u = zeros(Nx, Na);
center_x = Lx / 2;
center_a = La / 2;
sigma_x = Lx / 10;
sigma_a = La / 10;

u(1,1) = 1;

C = @(xx, tt) xx/Lx;

for ix = 1:Nx
    ix_plus = ix + 1;
    ix_minus = ix - 1;
    if ix == 1;
        ix_minus = ix;
    elseif ix == Nx
        ix_plus = ix ;
    end
    dc1_dx(ix) = (C(ix_plus*Dx,1) - C(ix*Dx,1) )/Dx;
    dc2_dx(ix) = (C(ix*Dx,1) - C(ix_minus*Dx,1))/Dx;
end



figure;
% Time stepping loop
for n = 1:Nt
    
    u_new = u;
    t_current = t(n);

    for ix = 1:Nx
        ix_plus = ix + 1;
        ix_minus = ix - 1;
        if ix == 1;
            ix_minus = ix;
        elseif ix == Nx
            ix_plus = ix ;
        end

        for ia = 1:Na
            ia_plus = ia + 1;
            ia_minus = ia - 1;
            if ia == 1;
                ia_minus = ia;
            elseif ia == Na
                ia_plus = ia ;
            end
            u_current = u(ix, ia);
            a_current = a(ia);

            % Diffusion term 
            x_laplacian = 0;
            x_laplacian = (u(ix_plus, ia) - 2 * u_current + u(ix_minus, ia)) / (Dx^2);

            % Chemotaxis term
            x_divergence = 0;
            u_plus_x = (u(ix_plus, ia) + u_current) / 2;
            u_minus_x = (u(ix_minus, ia) + u_current) / 2;

            x_divergence = (u_plus_x *  dc1_dx(ix) - u_minus_x * dc2_dx(ix)) / Dx;


            % 'a' diffusion term: Da_diffusion * d^2u/da^2 with reflective BCs
            a_laplacian = 0;
            a_laplacian = (u(ix, ia_plus) - 2 * u_current + u(ix, ia_minus )) / (Da^2);

            % 'a' derivative term: - d/da (f(u,a) * u)
            % a_derivative = 0;
            % 
            % 
            % if ia == Na % Absorbing boundary at ia = Na
            %     g_minus = f(u(ix, ia - 1), a(ia - 1)) * u(ix, ia - 1);
            %     g = f(u_current, a_current) * u_current;
            %     a_derivative = (g - g_minus) / Da;
            % 
            % elseif ia == 1 % Reflective boundary at ia = 1 (zero gradient of g)
            %     g = f(u_current, a_current) * u_current;
            %     g_plus = f(u(ix, ia + 1), a(ia + 1)) * u(ix, ia + 1);
            %     a_derivative = (g_plus - g) / Da;
            % 
            % else % Interior points
            %     g_plus = f(u(ix, ia + 1), a(ia + 1)) * u(ix, ia + 1);
            %     g_minus = f(u(ix, ia - 1), a(ia - 1)) * u(ix, ia - 1);
            %     a_derivative = (g_plus - g_minus) / (2 * Da);
            % end

            % 'a' derivative term: - (u * df/da + f * du/da)
            a_derivative = 0;
            u_current = u(ix, ia);
            a_current = a(ia);

            if ia > 1 && ia < Na % Interior points
                u_plus_a = u(ix, ia + 1);
                u_minus_a = u(ix, ia - 1);
                a_plus = a(ia + 1);
                a_minus = a(ia - 1);

                % Numerical approximation of df/da
                dfda = (f(u_plus_a, a_plus) - f(u_minus_a, a_minus)) / (2 * Da);

                % Numerical approximation of du/da
                duda = (u_plus_a - u_minus_a) / (2 * Da);

                a_derivative = -(u_current * dfda + f(u_current, a_current) * duda);

            elseif ia == Na % Absorbing boundary at ia = Na (using backward difference)
                u_minus_a = u(ix, ia - 1);
                a_minus = a(ia - 1);

                dfda = (f(u_current, a_current) - f(u_minus_a, a_minus)) / Da;
                duda = (u_current - u_minus_a) / Da;

                a_derivative = -(u_current * dfda + f(u_current, a_current) * duda);

            elseif ia == 1 % Reflective boundary at ia = 1 (using forward difference)
                u_plus_a = u(ix, ia + 1);
                a_plus = a(ia + 1);

                dfda = (f(u_plus_a, a_plus) - f(u_current, a_current)) / Da;
                duda = (u_plus_a - u_current) / Da;

                a_derivative = -(u_current * dfda + f(u_current, a_current) * duda);
                
            end


            % Update u using forward Euler
            u_new(ix, ia) = u_current + Dt * (D_x * x_laplacian - chi_x * x_divergence - chi_a * a_derivative + D_a * a_laplacian);
        end
    end

    u = u_new;

    if mod(n,10) == 0
        clf;
        plot(a , sum(u,1))


        % 
        % u_slice = squeeze(u(:, :));
        % surf(x, a, u_slice', 'EdgeColor','none');
        % % imagesc(x, a, u_slice');
        % 
        % colorbar;
        title(['Time = ', num2str(t_current, '%.2f'), ' mass = ', num2str(sum(u(:)), '%.2f')]);
        % xlabel('x');
        % ylabel('a');
        % colormap('jet');
        % view(2)

        drawnow
    end
end
