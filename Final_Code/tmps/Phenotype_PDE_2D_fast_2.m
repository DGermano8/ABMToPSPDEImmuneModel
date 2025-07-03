
plt_t = 1000;

% Define parameters
Np = 500;

Lx = 5; 
Ly = 5; 
La = 5; 
T  = 12;

Dx = 0.2; 
Dy = 0.2; 
Da = 0.1; 
Dt = 0.0001;

Nt = T/Dt;
Na = La/Da;
Nx = Lx/Dx;
Ny = Ly/Dy;

D = 0.10;
chi = 0.2; 
Da_diffusion = 0.01;
chia = 0.1; % 'a'-chemotaxis sensitivity
rho_plus = 0.5;
rho_minus = 0.1;

% Create grid
x = linspace(0, Lx, Nx); 
y = linspace(0, Ly, Ny); 
a = linspace(0, La, Na); 
t = linspace(0, T, Nt + 1);

% Initialize u (example initial condition)
u = zeros(Nx, Ny, Na);
% u(1,1,1) = 1;
% u(1,round(0.5*Ny),1) = 1;
% u(1,end,1) = 1;
u(1,:,round(0.1*Na)) = 1;
u = u/sum(u(:));

activation_proportion = zeros(Nt + 1,1);

% Initialize C (example function of x and t)
% C = @(xx, yy, tt) xx/Lx;
[Y,X, A] = meshgrid(y,x,a);
C = X/Lx;

% Define f(u, a) (example function)
Ic = zeros(Nx,Ny,Na);
Ic(round(0.5*Nx):end,:, :) = 1;

%%

figure;
% Time stepping loop
for n = 1:Nt
    
    u_new = u;
    t_current = Dt*n;
    
    % compute derivatives using finite volume scheme and no-flux boundary
    % conditions

    u_pad_x = cat(1, u(1,:,:), u, u(end,:,:));
    C_pad_x = cat(1, C(1,:,:), C, C(end,:,:));
    % Fluxes in x-direction
    u_avg_plus_x = (u_pad_x(2:end, :, :) + u_pad_x(1:(end-1),:,:)) / 2;
    dudx_plus = (u_pad_x(2:end, :,:) - u_pad_x(1:(end-1),:,:)) / Dx;
    dCdx_plus = (C_pad_x(2:end, :, :) - C_pad_x(1:(end-1),:, :)) / Dx;
    flux_x_plus = D .* dudx_plus - u_avg_plus_x .* chi .* dCdx_plus;

    u_avg_minus_x = (u_pad_x(2:end,:, :) + u_pad_x(1:(end-1),:, :)) / 2;
    dudx_minus = (u_pad_x(2:end,:, :) - u_pad_x(1:(end-1),:, :)) / Dx;
    dCdx_minus = (C_pad_x(2:end, :, :) - C_pad_x(1:(end-1),: ,:)) / Dx;
    flux_x_minus = D .* dudx_minus - u_avg_minus_x .* chi .* dCdx_minus;

    div_x_final = (flux_x_plus(2:end,:,:) - flux_x_minus(1:(end-1),:,:)) / Dx;

    
    u_pad_y = cat(2, u(:,1,:), u, u(:,end,:));
    C_pad_y = cat(2, C(:,1,:), C, C(:,end,:));
    % Fluxes in y-direction
    u_avg_plus_y = (u_pad_y(:, 2:end, :) + u_pad_y(:,1:(end-1),:)) / 2;
    dudy_plus = (u_pad_y(:, 2:end,:) - u_pad_y(:,1:(end-1),:)) / Dy;
    dCdy_plus = (C_pad_y(:, 2:end,:) - C_pad_y(:,1:(end-1),:)) / Dy;
    flux_y_plus = D .* dudy_plus - u_avg_plus_y .* chi .* dCdy_plus;

    u_avg_minus_y = (u_pad_y(:,2:end, :) + u_pad_y(:,1:(end-1), :)) / 2;
    dudy_minus = (u_pad_y(:,2:end, :) - u_pad_y(:,1:(end-1), :)) / Dy;
    dCdy_minus = (C_pad_y(:,2:end,:) - C_pad_y(:, 1:(end-1),:)) / Dy;
    flux_y_minus = D .* dudy_minus - u_avg_minus_y .* chi .* dCdy_minus;

    div_y_final = (flux_y_plus(:,2:end,:) - flux_y_minus(:,1:(end-1),:)) / Dy;    
    

    % u_pad_a = cat(3, u(:,:,2), u, u(:,:,end-1));
    % A_pad_a = cat(3, A(:,:,2), A, A(:,:,end-1));
    % I_pad_a = cat(3, Ic(:,:,1), Ic, Ic(:,:,end));
    % % flux is a-dirction, by upwinding
    % Flux = zeros(Nx,Ny,Na+2);
    % u_avg = u_pad_a(:,:,2:end-1);
    % a_i = A;
    % duda = (u_pad_a(:,:,2:end-1) - u_pad_a(:,:,1:(end-2))) / Da;
    % df = (rho_plus .* (1 - a_i./La) ) ;
    % Flux(:,:,2:end-1) = chia .* df .* u_avg - Da_diffusion .* duda;
    % flux_a_plus = Flux(:,:,3:end);
    % flux_a_minus = Flux(:,:,2:end-1);
    % 
    % div_a_final = (flux_a_plus - flux_a_minus)/Da;

    u_padded_a = zeros(Nx, Ny, Na + 2);
    u_padded_a(:,:,2:Na+1) = u; % Original u in the inner part

    a_padded = [a(1) a a(end)];
    
    % Update ghost cells for no-flux BC in a-direction
    f_a_0 = rho_plus * (1 - a(1) / La);
    f_a_La = rho_plus * (1 - a(end) / La);
    
    u_padded_a(:,:,1) = u_padded_a(:,:,3) - 2 * Da * chia * f_a_0 .* u_padded_a(:,:,2) / Da_diffusion;
    u_padded_a(:,:,Na+2) = u_padded_a(:,:,Na) + 2 * Da * chia * f_a_La .* u_padded_a(:,:,Na+1) / Da_diffusion;

    Flux_a = zeros(Nx, Ny, Na + 1);
    for j = 1:Na + 1
        % Location of the interface between cell j and j+1 (using padded indexing)
        a_interface = (a_padded(j) + a(min(j + 1, Na))) / 2; % Approximate interface location

        % Evaluate f at the interface (or use cell-centered value)
        df_interface = rho_plus * (1 - a_interface / La);
        
        % Upwind scheme for advection
        if chia * df_interface >= 0 % Advection in positive a direction
            u_upwind = u_padded_a(:,:,j+1); % Use value from the cell j (upstream)
        else % Advection in negative a direction (should not happen here if f >= 0)
            u_upwind = u_padded_a(:,:,j+2); % Use value from the cell j+1 (upstream)
        end
        
        % Flux at interface j
        Flux_a(:,:,j) = chia .* df_interface .* u_upwind - Da_diffusion .* (u_padded_a(:,:,j+2) - u_padded_a(:,:,j+1)) / Da;
    end
    
    div_a_final = (Flux_a(:,:,2:end) - Flux_a(:,:,1:end-1)) / Da;

    % Update u
    u = u + Dt * (div_x_final + div_y_final - div_a_final);


    activation_proportion(n) = squeeze(sum(u,[1,2]))'*a'/La;

    if true && mod(n,plt_t)==0
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
        axis([0 max(a) 0 max(squeeze(sum(u,[1,2])))])

    
        sgtitle(['time = ', num2str(Dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
    end
end
