
% Define parameters
Np = 500;

Lx = 5; 
Ly = 5; 
La = 10; 
T  = 100;

Dx = 0.2; 
Dy = 0.2; 
Da = 0.02; 
Dt = 0.005;

Nt = T/Dt;
Na = La/Da;
Nx = Lx/Dx;
Ny = Ly/Dy;

D = 0.5;
chi = 1; 
Da_diffusion = 0.001;
chia = 1 % 'a'-chemotaxis sensitivity
rho_plus = 0.9;
rho_minus = 0.2;

% Create grid
x = linspace(0, Lx, Nx); 
y = linspace(0, Ly, Ny); 
a = linspace(0, La, Na); 
t = linspace(0, T, Nt + 1);

% Initialize u (example initial condition)
u = zeros(Nx, Ny, Na);
u(1,:,1) = 1;
u = u/sum(u(:));

activation_proportion = zeros(Nt + 1,1);

% Initialize C (example function of x and t)
% C = @(xx, yy, tt) xx/Lx;
[Y,X, A] = meshgrid(y,x,a);
C = X/Lx;

% Define f(u, a) (example function)
Ic = zeros(Nx,Ny,Na);
Ic(round(0.5*Nx):end,:, :) = 1;


figure;
% Time stepping loop
for n = 1:Nt
    
    u_new = u;
    t_current = Dt*n;

    % for ia = 1:Na
    %     % Fluxes in x-direction
    %     u_avg_plus_x = (u(2:end, :, ia) + u(1:(end-1),:,ia)) / 2;
    %     dudx_plus = (u(2:end, :,ia) - u(1:(end-1),:,ia)) / Dx;
    %     dCdx_plus = (C(2:end, :) - C(1:(end-1),:)) / Dx;
    %     flux_x_plus = D .* dudx_plus - u_avg_plus_x .* chi .* dCdx_plus;
    %     flux_x_plus = [flux_x_plus;  zeros(1, Ny)];
    % 
    %     u_avg_minus_x = (u(2:end,:, ia) + u(1:(end-1),:, ia)) / 2;
    %     dudx_minus = (u(2:end,:, ia) - u(1:(end-1),:, ia)) / Dx;
    %     dCdx_minus = (C(2:end, :) - C(1:(end-1),:)) / Dx;
    %     flux_x_minus = D .* dudx_minus - u_avg_minus_x .* chi .* dCdx_minus;
    %     flux_x_minus = [zeros(1, Ny); flux_x_minus];
    % 
    %     div_x_final = (flux_x_plus - flux_x_minus) / Dx;
    % 
    %     % Fluxes in y-direction
    %     u_avg_plus_y = (u(:, 2:end, ia) + u(:,1:(end-1),ia)) / 2;
    %     dudy_plus = (u(:, 2:end,ia) - u(:,1:(end-1),ia)) / Dy;
    %     dCdy_plus = (C(:, 2:end) - C(:,1:(end-1))) / Dy;
    %     flux_y_plus = D .* dudy_plus - u_avg_plus_y .* chi .* dCdy_plus;
    %     flux_y_plus = [flux_y_plus, zeros(Nx, 1)];
    % 
    %     u_avg_minus_y = (u(:,2:end, ia) + u(:,1:(end-1), ia)) / 2;
    %     dudy_minus = (u(:,2:end, ia) - u(:,1:(end-1), ia)) / Dy;
    %     dCdy_minus = (C(:,2:end) - C(:, 1:(end-1))) / Dy;
    %     flux_y_minus = D .* dudy_minus - u_avg_minus_y .* chi .* dCdy_minus;
    %     flux_y_minus = [zeros(Nx, 1), flux_y_minus];
    % 
    %     div_y_final = (flux_y_plus - flux_y_minus) / Dy;
    % 
    % 
    %     % Fluxes in a-direction
    %     flux_a_plus = 0;
    %     flux_a_minus = 0;
    %     if ia < Na
    %         u_avg_plus_a = (u(:, :, ia + 1) + u(:, :, ia)) / 2;
    %         a_plus_half_a = a(ia) + Da / 2;
    %         duda_plus = (u(:, :, ia + 1) - u(:, :, ia)) / Da;
    %         flux_a_plus =  chia .* (rho_plus .* (1 - a_plus_half_a./La) .* Ic(:,:) - rho_minus .* a_plus_half_a./La .* (1-Ic(:,:)) ) .* u_avg_plus_a - Da_diffusion .* duda_plus;
    %     end
    %     if ia > 1
    %         u_avg_minus_a = (u(:, :, ia) + u(:, :, ia - 1)) / 2;
    %         a_minus_half_a = a(ia) - Da / 2;
    %         duda_minus = (u(:, :, ia) - u(:, :, ia - 1)) / Da;
    %         flux_a_minus = chia .* (rho_plus .* (1 - a_minus_half_a./La) .* Ic(:,:) - rho_minus .* a_minus_half_a./La .* (1-Ic(:,:)) ) .* u_avg_minus_a - Da_diffusion .* duda_minus;
    %     end
    % 
    %     % Boundary conditions for a
    %     if ia == 1
    %         % Outflow: flux at a=0 depends on the gradient into the domain
    %         % flux_a_minus = chia * f((u(ix, ia + 1) + u(ix, ia)) / 2, a(ia) + Da / 2) * (u(ix, ia + 1) + u(ix, ia)) / 2 - Da_diffusion * (u(ix, ia + 1) - u(ix, ia)) / Da;
    %         flux_a_minus = 0;
    %     elseif ia == Na
    %         % Absorbing: flux at a=La uses the value at the boundary
    %         % flux_a_plus = chia * f(u(ix, ia), a(ia) + Da / 2) * u(ix, ia) - Da_diffusion * (0 - u(ix, ia)) / Da; % Assume u outside is 0
    %         flux_a_plus = 0;
    %     end
    % 
    %     div_a_final = (flux_a_plus - flux_a_minus) / Da;
    % 
    %     % Update u
    %     u_new(:, :, ia) = u(:, :, ia) + Dt * (div_x_final + div_y_final - div_a_final);
    % end


    % Fluxes in x-direction
    u_avg_plus_x = (u(2:end, :, :) + u(1:(end-1),:,:)) / 2;
    dudx_plus = (u(2:end, :,:) - u(1:(end-1),:,:)) / Dx;
    dCdx_plus = (C(2:end, :, :) - C(1:(end-1),:, :)) / Dx;
    flux_x_plus = D .* dudx_plus - u_avg_plus_x .* chi .* dCdx_plus;
    flux_x_plus = [flux_x_plus;  zeros(1, Ny, Na)];

    u_avg_minus_x = (u(2:end,:, :) + u(1:(end-1),:, :)) / 2;
    dudx_minus = (u(2:end,:, :) - u(1:(end-1),:, :)) / Dx;
    dCdx_minus = (C(2:end, :, :) - C(1:(end-1),: ,:)) / Dx;
    flux_x_minus = D .* dudx_minus - u_avg_minus_x .* chi .* dCdx_minus;
    flux_x_minus = [zeros(1, Ny, Na); flux_x_minus];

    div_x_final = (flux_x_plus - flux_x_minus) / Dx;


    % Fluxes in y-direction
    u_avg_plus_y = (u(:, 2:end, :) + u(:,1:(end-1),:)) / 2;
    dudy_plus = (u(:, 2:end,:) - u(:,1:(end-1),:)) / Dy;
    dCdy_plus = (C(:, 2:end,:) - C(:,1:(end-1),:)) / Dy;
    flux_y_plus = D .* dudy_plus - u_avg_plus_y .* chi .* dCdy_plus;
    flux_y_plus = [flux_y_plus, zeros(Nx, 1, Na)];

    u_avg_minus_y = (u(:,2:end, :) + u(:,1:(end-1), :)) / 2;
    dudy_minus = (u(:,2:end, :) - u(:,1:(end-1), :)) / Dy;
    dCdy_minus = (C(:,2:end,:) - C(:, 1:(end-1),:)) / Dy;
    flux_y_minus = D .* dudy_minus - u_avg_minus_y .* chi .* dCdy_minus;
    flux_y_minus = [zeros(Nx, 1, Na), flux_y_minus];

    div_y_final = (flux_y_plus - flux_y_minus) / Dy;
    

    % Fluxes in a-direction
    u_avg_plus_a = (u(:, :, 2:end) + u(:, :, 1:(end-1))) / 2;
    a_plus_half_a = A(:,:, 1:(end-1)) + Da / 2;
    duda_plus = (u(:, :, 2:end) - u(:, :, 1:(end-1))) / Da;
    flux_a_plus =  chia .* (rho_plus .* (1 - a_plus_half_a./La) .* Ic(:,:, 1:(end-1)) - rho_minus .* a_plus_half_a./La .* (1-Ic(:,:, 1:(end-1))) ) .* u_avg_plus_a - Da_diffusion .* duda_plus;
    flux_a_plus = cat(3,flux_a_plus, zeros(Nx, Ny, 1));
    
    u_avg_minus_a = (u(:, :, 2:end) + u(:, :, 1:(end-1))) / 2;
    a_minus_half_a = A(:,:,2:end) - Da / 2;
    duda_minus = (u(:, :, 2:end) - u(:, :, 1:(end-1))) / Da;
    flux_a_minus = chia .* (rho_plus .* (1 - a_minus_half_a./La) .* Ic(:,:, 2:end) - rho_minus .* a_minus_half_a./La .* (1-Ic(:,:, 2:end)) ) .* u_avg_minus_a - Da_diffusion .* duda_minus;
    flux_a_minus = cat(3, zeros(Nx, Ny, 1),flux_a_minus);

    div_a_final = (flux_a_plus - flux_a_minus) / Da;
    
    % Update u
    u = u + Dt * (div_x_final + div_y_final - div_a_final);
    u = u/sum(u(:));

    activation_proportion(n) = squeeze(sum(u,[1,2]))'*a'/La;

    if true && mod(n,100)==0
        clf;
    
        subplot(2,2,1)
        u_slice = squeeze(sum(u,3));
        surf(x, y, u_slice', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        colormap('jet');
        view(2)
        axis([0 Lx 0 Ly])
    
        av_a = zeros(Nx, Ny);
        for ix=1:Nx
            for iy=1:Ny
                av_a(ix,iy) = a*squeeze(u(ix,iy,:));
            end
        end
        av_a = av_a/La*sum(u(:));

    
        subplot(2,2,2)
        surf(x, y, av_a', 'EdgeColor','none');
        colorbar;
        xlabel('x');
        ylabel('y');
        colormap('jet');
        view(2)
        axis([0 Lx 0 Ly])
        
        
        subplot(2,2,3)
        plot(1:n,activation_proportion(1:n))
        axis([0 Nt 0 1])

        subplot(2,2,4)
        plot(a, squeeze(sum(u,[1,2])))

    
        sgtitle(['time = ', num2str(Dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end
end
