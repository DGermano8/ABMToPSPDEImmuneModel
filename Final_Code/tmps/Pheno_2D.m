% Define parameters
Nx = 101; Ny = 101; Na = 30; Nt = 1000;
Lx = 10; Ly = 10; La = 5; T = 10;
Dx = Lx / (Nx - 1); Dy = Ly / (Ny - 1); Da = La / (Na - 1); Dt = T / Nt;

D = 0.2; 
chi = 0.1; 
Da_diffusion = 0.0;
chia = 0.00;

% Create grid
x = linspace(0, Lx, Nx); 
y = linspace(0, Ly, Ny); 
a = linspace(0, La, Na); 
t = linspace(0, T, Nt + 1);

% Initialize u (example initial condition)
u = zeros(Nx, Ny, Na);
center_x = Lx / 2; center_y = Ly / 2; center_a = La / 2;
sigma_x = Lx / 10; sigma_y = Ly / 10; sigma_a = La / 10;
% for ix = 1:Nx
%     for iy = 1:Ny
%         for ia = 1:Na
%             u(ix, iy, ia) = exp(-((x(ix) - center_x)^2) / (2 * sigma_x^2) - ((y(iy) - center_y)^2) / (2 * sigma_y^2) - ((a(ia) - center_a)^2) / (2 * sigma_a^2));
%         end
%     end
% end

u(1,50) = 1;

% Initialize C (example function of x, y, and t)
C = @(xx, yy, tt) xx/Lx;
dCdx = @(xx, yy, tt) 1/Lx;
dCdy = @(xx, yy, tt) 0;

% Define f(u, a) (example function)
f = @(uu, aa) 1 + 0.1 * uu .* aa;

% Time stepping loop
for n = 1:Nt
    u_new = u;
    t_current = t(n);

    for ix = 1:Nx
        for iy = 1:Ny
            for ia = 1:Na
                % Spatial Diffusion in x
                laplacian_u_x = 0;
                if ix == 1
                    % laplacian_u_x = (u(ix + 1, iy, ia) - 2 * u(ix, iy, ia) + u(ix + 1, iy, ia)) / (Dx^2); % Reflective
                    laplacian_u_x = (u(ix, iy, ia) - 2 * u(ix, iy, ia) + u(ix + 1, iy, ia)) / (Dx^2); % Reflective
                elseif ix == Nx
                    % laplacian_u_x = (u(ix - 1, iy, ia) - 2 * u(ix, iy, ia) + u(ix - 1, iy, ia)) / (Dx^2); % Reflective
                    laplacian_u_x = (u(ix - 1, iy, ia) - 2 * u(ix, iy, ia) + u(ix, iy, ia)) / (Dx^2); % Reflective
                else
                    laplacian_u_x = (u(ix + 1, iy, ia) - 2 * u(ix, iy, ia) + u(ix - 1, iy, ia)) / (Dx^2);
                end

                % Spatial Diffusion in y
                laplacian_u_y = 0;
                if iy == 1
                    % laplacian_u_y = (u(ix, iy + 1, ia) - 2 * u(ix, iy, ia) + u(ix, iy + 1, ia)) / (Dy^2); % Reflective
                    laplacian_u_y = (u(ix, iy, ia) - 2 * u(ix, iy, ia) + u(ix, iy + 1, ia)) / (Dy^2); % Reflective
                elseif iy == Ny
                    % laplacian_u_y = (u(ix, iy - 1, ia) - 2 * u(ix, iy, ia) + u(ix, iy - 1, ia)) / (Dy^2); % Reflective
                    laplacian_u_y = (u(ix, iy - 1, ia) - 2 * u(ix, iy, ia) + u(ix, iy, ia)) / (Dy^2); % Reflective
                else
                    laplacian_u_y = (u(ix, iy + 1, ia) - 2 * u(ix, iy, ia) + u(ix, iy - 1, ia)) / (Dy^2);
                end

                % Chemotaxis in x
                chemotaxis_x = 0;
                if ix == 1
                    u_plus_x = (u(ix + 1, iy, ia) + u(ix, iy, ia)) / 2;
                    u_minus_x = (u(ix, iy, ia) + u(ix, iy, ia)) / 2;
                    chemotaxis_x = (u_plus_x  * dCdx(x(ix), y(iy), t_current) - u_minus_x  * dCdx(x(ix), y(iy), t_current)) / Dx;
                    % u_plus_x = (u(ix, iy, ia) + u(ix, iy, ia)) / 2;
                    % chemotaxis_x = (u_plus_x *  dCdx(x(ix), y(iy), t_current) - u_plus_x *  dCdx(x(ix), y(iy), t_current)) / Dx;
                elseif ix == Nx
                    u_plus_x = (u(ix, iy, ia) + u(ix, iy, ia)) / 2;
                    u_minus_x = (u(ix - 1, iy, ia) + u(ix, iy, ia)) / 2;
                    chemotaxis_x = (u_plus_x  * dCdx(x(ix), y(iy), t_current) - u_minus_x *  dCdx(x(ix), y(iy), t_current)) / Dx;
                    % u_minus_x = (u(ix, iy, ia) + u(ix, iy, ia)) / 2;
                    % chemotaxis_x = (u_minus_x *  dCdx(x(ix), y(iy), t_current) - u_minus_x *  dCdx(x(ix), y(iy), t_current)) / Dx;
                else
                    u_plus_x = (u(ix + 1, iy, ia) + u(ix, iy, ia)) / 2;
                    u_minus_x = (u(ix, iy, ia) + u(ix - 1, iy, ia)) / 2;
                    chemotaxis_x =  (u_plus_x  * dCdx(x(ix + 1), y(iy), t_current) - u_minus_x *  dCdx(x(ix - 1), y(iy), t_current)) / Dx;
                end
                chemotaxis_x = -chi *chemotaxis_x; % Sign correction for the divergence term

                % Chemotaxis in y
                chemotaxis_y = 0;
                if iy == 1
                    % u_plus_y = (u(ix, iy + 1, ia) + u(ix, iy, ia)) / 2;
                    % chemotaxis_y = (u_plus_y * dCdy(x(ix), y(iy + 1), t_current) - u_plus_y * dCdy(x(ix), y(iy), t_current)) / Dy;
                    u_plus_y = (u(ix, iy, ia) + u(ix, iy, ia)) / 2;
                    chemotaxis_y = (u_plus_y * dCdy(x(ix), y(iy), t_current) - u_plus_y *  dCdy(x(ix), y(iy), t_current)) / Dy;
                elseif iy == Ny
                    % u_minus_y = (u(ix, iy - 1, ia) + u(ix, iy, ia)) / 2;
                    % chemotaxis_y = (u_minus_y * dCdy(x(ix), y(iy), t_current) - u_minus_y *  dCdy(x(ix), y(iy - 1), t_current)) / Dy;
                    u_minus_y = (u(ix, iy, ia) + u(ix, iy, ia)) / 2;
                    chemotaxis_y = (u_minus_y * dCdy(x(ix), y(iy), t_current) - u_minus_y *  dCdy(x(ix), y(iy), t_current)) / Dy;
                else
                    u_plus_y = (u(ix, iy + 1, ia) + u(ix, iy, ia)) / 2;
                    u_minus_y = (u(ix, iy, ia) + u(ix, iy - 1, ia)) / 2;
                    chemotaxis_y = (u_plus_y  * dCdy(x(ix), y(iy + 1), t_current) - u_minus_y *  dCdy(x(ix), y(iy - 1), t_current)) / Dy;
                end
                chemotaxis_y = -chi *chemotaxis_y; % Sign correction for the divergence term

                % 'a' Derivative and Diffusion
                div_a = 0;
                if ia > 1 && ia < Na
                    u_plus_a = (u(ix, iy, ia + 1) + u(ix, iy, ia)) / 2;
                    F_plus_a = chia * f(u_plus_a, a(ia) + Da / 2) * u_plus_a - Da_diffusion * (u(ix, iy, ia + 1) - u(ix, iy, ia)) / Da;
                    u_minus_a = (u(ix, iy, ia) + u(ix, iy, ia - 1)) / 2;
                    F_minus_a = chia * f(u_minus_a, a(ia) - Da / 2) * u_minus_a - Da_diffusion * (u(ix, iy, ia) - u(ix, iy, ia - 1)) / Da;
                    div_a = (F_plus_a - F_minus_a) / Da;
                elseif ia == 1 % Outflow
                    u_plus_a = (u(ix, iy, ia + 1) + u(ix, iy, ia)) / 2;
                    F_plus_a = chia * f(u_plus_a, a(ia) + Da / 2) * u_plus_a - Da_diffusion * (u(ix, iy, ia + 1) - u(ix, iy, ia)) / Da;
                    % Approximate flux at ia=1/2 based on values at ia=1 and ia=2
                    u_boundary = u(ix, iy, ia);
                    F_minus_a = chia * f(u_boundary, a(ia) - Da / 2) * u_boundary - Da_diffusion * (u_boundary - u(ix, iy, ia + 1)) / Da; % Simplified
                    div_a = (F_plus_a - F_minus_a) / Da;
                elseif ia == Na % Absorbing
                    u_minus_a = (u(ix, iy, ia) + u(ix, iy, ia - 1)) / 2;
                    F_minus_a = chia * f(u_minus_a, a(ia) - Da / 2) * u_minus_a - Da_diffusion * (u(ix, iy, ia) - u(ix, iy, ia - 1)) / Da;
                    % Approximate flux at ia=Na+1/2 based on value at ia=Na
                    u_boundary = u(ix, iy, ia);
                    F_plus_a = chia * f(u_boundary, a(ia) + Da / 2) * u_boundary - Da_diffusion * (0 - u_boundary) / Da; % Absorbing u=0 outside
                    div_a = (F_plus_a - F_minus_a) / Da;
                end
                div_a = 0;

                % Update u
                u_new(ix, iy, ia) = u(ix, iy, ia) + Dt * (D * (laplacian_u_x + laplacian_u_y) - (chemotaxis_x + chemotaxis_y) - div_a);
            end
        end
    end
    u = u_new;

if true && mod(n,1)==0
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
    surf(x, y, av_a, 'EdgeColor','none');
    colorbar;
    xlabel('x');
    ylabel('y');
    colormap('jet');
    view(2)
    


    sgtitle(['time = ', num2str(Dt*n), ' mass = ', num2str(sum(u(:)))])
    drawnow;
end



end


