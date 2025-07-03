
%%

% --- Example Usage in a separate script or command window ---
% Define parameters
L = 1.0;
T = 20.0;
Nx = 100;
D = 0.1;

% Initial condition: Gaussian pulse centered on the left
initial_gaussian = @(x) exp(-((x - 0.2 * L) / (0.1 * L)).^2);

% Run the simulation
[final_concentration, wall_conc_history] = diffusion_reflective_stickywall_matlab(L, T, Nx, D, initial_gaussian);

% Print final wall concentration
fprintf('Final concentration on the wall: %f\n', wall_conc_history(end));

% Optional: Plot wall concentration over time
time_points = linspace(0, T, length(wall_conc_history));
figure;
plot(time_points, wall_conc_history);
xlabel('Time');
ylabel('Concentration on Sticky Wall');
title('Accumulation on Sticky Wall over Time (Matlab)');



%%

function [u, wall_concentration] = diffusion_reflective_stickywall_matlab(L, T, Nx, D, initial_condition_func, plot_steps)
% Solves 1D diffusion with reflective (x=0) and sticky wall (x=L) BCs in Matlab.

    if nargin < 6
        plot_steps = 50; % Default plot steps if not provided
    end

    dx = L / (Nx - 1);
    dt = 0.5 * dx^2 / D;
    r = D * dt / dx^2;
    Nt = floor(T / dt) + 1; % Use floor to ensure integer number of steps

    x = linspace(0, L, Nx);
    u = initial_condition_func(x);
    disp(['initial mass: ', num2str(sum(u)*dx)])
    wall_concentration = zeros(1, Nt); % Array to store wall concentration over time

    % Initialize plot (optional, for visualization)
    figure;
    plot_line = plot(x, u);
    hold on;
    wall_line = plot(L, 0, 'ro'); % Plot wall concentration as a red dot on the right edge
    hold off;
    xlabel('x');
    ylabel('Concentration u(x, t)');
    title('Diffusion with Reflective & Sticky Wall (Matlab)');
    xlim([0, L*1.1]); % Extend x-axis to see wall point
    ylim([0, max(u) * 1.1]);
    drawnow; % Force initial plot to appear

    for j = 1:Nt % Matlab uses 1-based indexing
        u_new = u; % Initialize u_new with current u
        flux_to_wall = -D * (u(Nx) - u(Nx-1)) / dx; % Flux towards the wall (using current u)
        if j > 1
            wall_concentration(j) = wall_concentration(j-1) + dt * flux_to_wall; % Accumulate on wall
        else
            wall_concentration(j) = 0;
        end


        for i = 1:Nx
            if i == 1 % Reflective boundary at x=0 (Matlab index 1)
                u_new(i) = u(i) + 2 * r * (u(i+1) - u(i));
            elseif i == Nx % Sticky Wall boundary condition for *diffusion* calculation (Matlab index Nx)
                u_new(i) = 0.0; % Set concentration at boundary to 0 for next time step's diffusion
            elseif i == 2
                u_new(i) = u(i) + 2 * r * (u(i+1) - u(i));
            else % Interior points
                u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
            end
        end
        u = u_new;

        if mod(j, plot_steps) == 0
            set(plot_line, 'YData', u);
            set(wall_line, 'XData', L, 'YData', wall_concentration(j)); % Update wall point
            ylim([0, max([max(u) * 1.1, wall_concentration(j) * 1.1, 0.1])]); % Adjust y-limit, ensure non-zero lower bound for ylim
            title(['Diffusion with Reflective & Sticky Wall (Matlab), time: ', num2str(j*dt), ' mass: ', num2str(sum(u))]);
            drawnow; % Update plot
        end
    end

    hold off; % Release hold on plot
end
