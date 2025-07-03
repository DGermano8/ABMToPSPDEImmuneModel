clear; clc;

%% Parameters
Nx = 50; Na = 20; % Grid points in x, y, and phenotype space
dx = 0.1; da = 0.05; % Grid spacing
D = 0.1; % Diffusion coefficient
chi = 0.5; % Taxis sensitivity
dt = 0.01; % Time step
T = 100; % Total simulation time

%% Initialize variables
u = zeros(Nx, Na); % T-cell density
C = ones(Nx,1); % Chemokine concentration (for taxis, modify as needed)
f = @(u, a) 0.1 * u.* (1 - a); % Activation function

u(1,1) = 1;
u = u/sum(u(:));

t = 0;
figure;

%%
for t  = 0:dt:T
    sum(u(:))
    u_new = u;
    
    % Diffusion (Central Difference)
    for i = 1:Nx
        for k = 1:Na
            if i > 1 && i < Nx
                dCdx = (C(i+1) - C(i-1)) / (2*dx);
                u_new(i,k) = u(i,k) + dt *( D * ((u(i+1,k) - 2*u(i,k) + u(i-1,k))/dx^2 ) - chi * (u(i,k) * dCdx)/dx );

            elseif i == 1

                dCdx = (C(i+1) - C(i)) / (2*dx);
                u_new(i,k) = u(i,k) + dt *( D * ((u(i+1,k) - 2*u(i,k) + u(i,k))/dx^2 ) - chi * (u(i,k) * dCdx)/dx );

            elseif i == Nx

                dCdx = (C(i) - C(i-1)) / (2*dx);
                u_new(i,k) = u(i,k) + dt *( D * ((u(i,k) - 2*u(i,k) + u(i-1,k))/dx^2 ) - chi * (u(i,k) * dCdx)/dx );

            end
        end
    end
    
    % Taxis (Upwind scheme)
    % for i = 2:Nx-1
    %     for k = 1:Na
    %         dCdx = (C(i+1) - C(i-1)) / (2*dx);
    % 
    %         u_new(i,k) = u_new(i,k) - dt * chi * (u(i,k) * dCdx)/dx;
    %     end
    % end
    
    % Advection in phenotype space (Upwind)
    for i = 1:Nx
        for k = 1:Na % Upwind for stability
            if k > 1
                u_new(i,k) = u_new(i,k) - dt * ( f(u(i,k), k*da) * (u(i,k) - u(i,k-1))/da);
            elseif k == 1
                u_new(i,k) = u_new(i,k) - dt * (f(u(i,k), k*da) * (u(i,k) - u(i,k+1))/da);
            end
        end
    end
    
    % Apply reflective BCs (Neumann conditions)
    % u_new(1,:) = u_new(2,:); u_new(end,:) = u_new(end-1,:);
    
    % % Reflective boundary conditions in x-direction
    % for k = 1:Na
    %     u_new(1,k) = u_new(2,k); 
    %     u_new(Nx,k) = u_new(Nx-1,k);
    % end
    
    if ismember(round(t,2), [0:1:T])
        % Visualization
        clf;
        subplot(1,3,1)
        plot(1:Nx, sum(u,2)'); % Summing over phenotype space
        title('T-cell Density Distribution');
        xlabel('Spatial position (x)');
        axis([0 Nx+1 0 1.1*max(sum(u,2)') ])

        subplot(1,3,2)
        plot(1:Na, sum(u,1)'); % Summing over phenotype space
        title('T-cell Density Distribution');
        xlabel('Spatial position (x)');
        axis([0 Na+1 0 1.1*max(sum(u,1)') ])

        subplot(1,3,3)
        [X, A] = meshgrid(1:Nx, (1:Na)*da);
        surf(X, A, u' ); % Summing over y-dimension
        xlabel('Spatial position (x)');
        ylabel('Phenotype (a)');
        zlabel('Density');
        title('Phenotype Density Distribution');
        shading interp;
        
        sgtitle(['Time = ', num2str(t)])
        drawnow;
    end

    u = u_new;
end

