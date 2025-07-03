clear; clc;

%% Parameters
Nx = 50; Ny = 60; Na = 20; % Grid points in x, y, and phenotype space
dx = 0.1; dy = 0.1; da = 0.05; % Grid spacing
D = 0.1; % Diffusion coefficient
chi = 0.5; % Taxis sensitivity
dt = 0.01; % Time step
T = 100; % Total simulation time

%% Initialize variables
u = zeros(Nx, Ny, Na); % T-cell density
C = ones(Nx, Ny); % Chemokine concentration (for taxis, modify as needed)
f = @(u, a) 0.1 * u .* (1 - a); % Activation function

u(1,:,:) = 1;
u = u/sum(u(:));

t = 0;
figure;

while t < T
    u_new = u;
    
    % Diffusion (Central Difference)
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 1:Na
                u_new(i,j,k) = u(i,j,k) + dt * D * ((u(i+1,j,k) - 2*u(i,j,k) + u(i-1,j,k))/dx^2 + ...
                                            (u(i,j+1,k) - 2*u(i,j,k) + u(i,j-1,k))/dy^2);
            end
        end
    end
    
    % Taxis (Upwind scheme)
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 1:Na
                dCdx = (C(i+1,j) - C(i-1,j)) / (2*dx);
                dCdy = (C(i,j+1) - C(i,j-1)) / (2*dy);
                
                u_new(i,j,k) = u_new(i,j,k) - dt * chi * ((u(i,j,k) * dCdx)/dx + (u(i,j,k) * dCdy)/dy);
            end
        end
    end
    
    % Advection in phenotype space (Upwind)
    for i = 1:Nx
        for j = 1:Ny
            for k = 2:Na % Upwind for stability
                u_new(i,j,k) = u_new(i,j,k) - dt * (f(u(i,j,k), k*da) * (u(i,j,k) - u(i,j,k-1))/da);
            end
        end
    end
    
    % Apply reflective BCs (Neumann conditions)
    u_new(1,:,:) = u_new(2,:,:); u_new(end,:,:) = u_new(end-1,:,:);
    u_new(:,1,:) = u_new(:,2,:); u_new(:,end,:) = u_new(:,end-1,:);
    
    % Reflective boundary conditions in x-direction
    for j = 1:Ny
        for k = 1:Na
            u_new(1,j,k) = u_new(2,j,k); 
            u_new(Nx,j,k) = u_new(Nx-1,j,k);
        end
    end
    
    if ismember(round(t,2), [0:5:T])
        % Visualization
        clf;
        subplot(1,2,1)
        [X, Y] = meshgrid(1:Nx, 1:Ny);
        surf(X, Y, sum(u,3)'); % Summing over phenotype space
        title('T-cell Density Distribution');
        xlabel('Spatial position (x)');
        ylabel('Spatial position (y)');
        shading interp;

        subplot(1,2,2)
        [X, A] = meshgrid(1:Nx, (1:Na)*da);
        surf(X, A, squeeze(sum(u(:,:, :), 2))' ); % Summing over y-dimension
        xlabel('Spatial position (x)');
        ylabel('Phenotype (a)');
        zlabel('Density');
        title('Phenotype Density Distribution');
        shading interp;
        
        sgtitle(['Time = ', num2str(t)])
        drawnow;
    end

    u = u_new;
    sum(u_new,'all')
    t = t + dt;
end

%%
figure;
[X, A] = meshgrid(1:Nx, (1:Na)*da);
surf(X, A, squeeze(sum(u(:,:, :), 2))'); % Summing over y-dimension
xlabel('Spatial position (x)');
ylabel('Phenotype (a)');
zlabel('Density');
title('3D Phenotype Density Distribution');
shading interp;


