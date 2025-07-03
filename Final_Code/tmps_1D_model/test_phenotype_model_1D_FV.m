% \frac{\partial u}{\partial t} = 
%       + \nabla   ( D \nabla u       - u \chi \nabla C )
%       - \nabla_a ( \chi_a f(u,a) u  - D_a  \nabla_a u )

    % Define parameters
    Nx = 200;          % Number of grid points in x
    Na = 100;           % Number of grid points in a
    Nt = 10000;         % Number of time steps
    Lx = 10;           % Length of spatial domain
    La = 5;            % Length of 'a' domain
    T = 10;            % Total time
    Dx = Lx / (Nx - 1); % Spatial step size
    Da = La / (Na - 1); % 'a' step size
    Dt = T / Nt;       % Time step size
    D = 0.01;           % Diffusion coefficient in x
    chi = 1;        % Chemotaxis sensitivity
    Da_diffusion = 0.01; % Diffusion coefficient in a
    chi_a = 0.0;

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
    % u = exp(-((x' - center_x).^2) / (2 * sigma_x^2)) .* exp(-((a - center_a).^2) / (2 * sigma_a^2));
    u(1,1) = 1;

    % Initialize C (example function of x and t)
    % C = @(xx, tt) 0.5 * (1 + sin(2 * pi * xx / Lx)) * exp(-tt / 5);
    % dCdx = @(xx, tt) 0.5 * (2 * pi / Lx) * cos(2 * pi * xx / Lx) * exp(-tt / 5);
    
    C = @(xx, tt) xx/Lx;
    dCdx = @(xx, tt) 1/Lx;
    
    % Define f(u, a) (example function)
    % f = @(uu, aa) 1 + 0.1 * uu .* aa;
    f = @(uu,aa) 1 * (1 - aa/La);

    % Store solution over time
    u_history = zeros(Nt + 1, Nx, Na);
    u_history(1, :, :) = u;

    % Time stepping loop
    for n = 1:Nt
        u_new = u; % Initialize for the next time step
        t_current = t(n);

        for ix = 1:Nx
            for ia = 1:Na
                % Spatial Diffusion (Finite Volume)
                laplacian_u_x = 0;
                if ix == 1
                    laplacian_u_x = (u(ix, ia) - 2 * u(ix, ia) + u(ix + 1, ia)) / (Dx^2); % Reflective
                elseif ix == Nx
                    laplacian_u_x = (u(ix, ia) - 2 * u(ix, ia) + u(ix - 1, ia)) / (Dx^2); % Reflective
                else
                    laplacian_u_x = (u(ix + 1, ia) - 2 * u(ix, ia) + u(ix - 1, ia)) / (Dx^2);
                end

                % Chemotaxis (Finite Volume)
                divergence_term = 0;
                if ix == 1
                    flux_plus = ((u(ix + 1, ia) + u(ix, ia)) / 2) * chi * dCdx(x(ix + 1), t_current);
                    flux_minus = ((u(ix + 1, ia) + u(ix, ia)) / 2) * chi * dCdx(x(ix), t_current); % Reflective
                    divergence_term = (flux_plus - flux_minus) / Dx;
                elseif ix == Nx
                    flux_plus = ((u(ix - 1, ia) + u(ix, ia)) / 2) * chi * dCdx(x(ix), t_current); % Reflective
                    flux_minus = ((u(ix - 1, ia) + u(ix, ia)) / 2) * chi * dCdx(x(ix - 1), t_current);
                    divergence_term = (flux_plus - flux_minus) / Dx;
                else
                    flux_plus = ((u(ix + 1, ia) + u(ix, ia)) / 2) * chi * dCdx(x(ix + 1), t_current);
                    flux_minus = ((u(ix, ia) + u(ix - 1, ia)) / 2) * chi * dCdx(x(ix - 1), t_current);
                    divergence_term = (flux_plus - flux_minus) / Dx;
                end

                % 'a' Derivative (Finite Volume)
                a_derivative_term = 0;
                if ia > 1 && ia < Na
                    u_plus_half_a = (u(ix, ia) + u(ix, ia + 1)) / 2;
                    F_plus_half_a = f(u_plus_half_a, a(ia) + Da / 2) * u_plus_half_a;

                    u_minus_half_a = (u(ix, ia) + u(ix, ia - 1)) / 2;
                    F_minus_half_a = f(u_minus_half_a, a(ia) - Da / 2) * u_minus_half_a;

                    a_derivative_term = (F_plus_half_a - F_minus_half_a) / Da; % Note the sign change here due to the -d/da term
                elseif ia == 1 % Reflective BC: F_{1/2} = 0
                    u_plus_half_a = (u(ix, ia) + u(ix, ia + 1)) / 2;
                    F_plus_half_a = f(u_plus_half_a, a(ia) + Da / 2) * u_plus_half_a;
                    a_derivative_term = (F_plus_half_a - 0) / Da; % Note the sign change
                elseif ia == Na % Absorbing BC: F_{Na+1/2} based on value at Na
                    u_minus_half_a = (u(ix, ia) + u(ix, ia - 1)) / 2;
                    F_minus_half_a = f(u_minus_half_a, a(ia) - Da / 2) * u_minus_half_a;
                    F_plus_half_a = f(u(ix, ia), a(ia) + Da / 2) * u(ix, ia);
                    a_derivative_term = (F_plus_half_a - F_minus_half_a) / Da; % Note the sign change
                end

                % 'a' Diffusion (Finite Volume)
                laplacian_u_a = 0;
                if ia == 1
                    laplacian_u_a = (u(ix, ia) - 2 * u(ix, ia) + u(ix, ia + 1)) / (Da^2); % Reflective
                elseif ia == Na
                    laplacian_u_a = (u(ix, ia) - 2 * u(ix, ia) + u(ix, ia - 1)) / (Da^2); % Reflective
                else
                    laplacian_u_a = (u(ix, ia + 1) - 2 * u(ix, ia) + u(ix, ia - 1)) / (Da^2);
                end

                % Update u
                u_new(ix, ia) = u(ix, ia) + Dt * (D * laplacian_u_x - divergence_term - chi_a* a_derivative_term + Da_diffusion * laplacian_u_a);
            end
        end
        u = u_new;
        % u_history(n + 1, :, :) =  u;

        if mod(n,10) == 0
            clf;
            % plot(a , sum(u,1))
    
    

            u_slice = squeeze(u(:, :));
            surf(x, a, u_slice', 'EdgeColor','none');
            % imagesc(x, a, u_slice');

            colorbar;
            title(['Time = ', num2str(t_current, '%.2f'), ' mass = ', num2str(sum(u(:)), '%.2f')]);
            xlabel('x');
            ylabel('a');
            colormap('jet');
            view(2)
    
            drawnow
        end

    end
    % 
    % % Visualize the solution
    % figure;
    % for i_t = 1:10:Nt + 1
    %     subplot(2, 5, ceil(i_t / 10));
    %     u_slice = squeeze(u_history(i_t, :, :));
    %     imagesc(x, a, u_slice');
    %     colorbar;
    %     title(['Time = ', num2str(t(i_t), '%.2f')]);
    %     xlabel('x');
    %     ylabel('a');
    %     colormap('jet');
    % end
    % sgtitle('Evolution of u(x, a, t) using Finite Volume Method with BCs');
