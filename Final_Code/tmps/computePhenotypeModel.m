function u = computePhenotypeModel(u, C, Ic, A_Ic, Ic_vals, params, A)
    

    Nx = params.Nx;
    Ny = params.Ny;
    Na = params.Na;

    La = params.La;


    Dx = params.dx;
    Dy = params.dy;
    Da = params.da;
    D = params.D;
    chi = params.C_chi;
    
    rho_plus = params.P_A;
    rho_minus = params.P_D;

    chia = 0;
    Da_diffusion = 0.0;
    

    u_avg_plus_x = (u(2:end, :, :) + u(1:(end-1),:,:)) / 2;
    dudx_plus = (u(2:end, :,:) - u(1:(end-1),:,:)) / Dx;
    dCdx_plus = (C(2:end, :, :) - C(1:(end-1),:, :)) / Dx;
    flux_x_plus = D .* dudx_plus - u_avg_plus_x .* chi .* dCdx_plus;
    flux_x_plus = [flux_x_plus;  zeros(1, Ny, Na)];
    flux_x_plus(Ic == 1) = 0;

    u_avg_minus_x = (u(2:end,:, :) + u(1:(end-1),:, :)) / 2;
    dudx_minus = (u(2:end,:, :) - u(1:(end-1),:, :)) / Dx;
    dCdx_minus = (C(2:end, :, :) - C(1:(end-1),: ,:)) / Dx;
    flux_x_minus = D .* dudx_minus - u_avg_minus_x .* chi .* dCdx_minus;
    flux_x_minus = [zeros(1, Ny, Na); flux_x_minus];
    flux_x_minus(Ic == 1) = 0;

    div_x_final = (flux_x_plus - flux_x_minus) / Dx;


    % Fluxes in y-direction
    u_avg_plus_y = (u(:, 2:end, :) + u(:,1:(end-1),:)) / 2;
    dudy_plus = (u(:, 2:end,:) - u(:,1:(end-1),:)) / Dy;
    dCdy_plus = (C(:, 2:end,:) - C(:,1:(end-1),:)) / Dy;
    flux_y_plus = D .* dudy_plus - u_avg_plus_y .* chi .* dCdy_plus;
    flux_y_plus = [flux_y_plus, zeros(Nx, 1, Na)];
    flux_y_plus(Ic == 1) = 0;

    u_avg_minus_y = (u(:,2:end, :) + u(:,1:(end-1), :)) / 2;
    dudy_minus = (u(:,2:end, :) - u(:,1:(end-1), :)) / Dy;
    dCdy_minus = (C(:,2:end,:) - C(:, 1:(end-1),:)) / Dy;
    flux_y_minus = D .* dudy_minus - u_avg_minus_y .* chi .* dCdy_minus;
    flux_y_minus = [zeros(Nx, 1, Na), flux_y_minus];
    flux_y_minus(Ic == 1) = 0;

    div_y_final = (flux_y_plus - flux_y_minus) / Dy;
    

    % Fluxes in a-direction
    u_avg_plus_a = (u(:, :, 2:end) + u(:, :, 1:(end-1))) / 2;
    a_plus_half_a = A(:,:, 1:(end-1)) + Da / 2;
    duda_plus = (u(:, :, 2:end) - u(:, :, 1:(end-1))) / Da;
    flux_a_plus =  chia .* (rho_plus .* (1 - a_plus_half_a./La) .* A_Ic(:,:, 1:(end-1)) - rho_minus .* a_plus_half_a./La .* (1-A_Ic(:,:, 1:(end-1))) ) .* u_avg_plus_a - Da_diffusion .* duda_plus;
    flux_a_plus = cat(3,flux_a_plus, zeros(Nx, Ny, 1));
    
    u_avg_minus_a = (u(:, :, 2:end) + u(:, :, 1:(end-1))) / 2;
    a_minus_half_a = A(:,:,2:end) - Da / 2;
    duda_minus = (u(:, :, 2:end) - u(:, :, 1:(end-1))) / Da;
    flux_a_minus = chia .* (rho_plus .* (1 - a_minus_half_a./La) .* A_Ic(:,:, 2:end) - rho_minus .* a_minus_half_a./La .* (1-A_Ic(:,:, 2:end)) ) .* u_avg_minus_a - Da_diffusion .* duda_minus;
    flux_a_minus = cat(3, zeros(Nx, Ny, 1),flux_a_minus);

    div_a_final = (flux_a_plus - flux_a_minus) / Da;

    div_a_final = zeros(Nx,Ny,Na);
    
    % Update u
    u = u + params.dt * (div_x_final + div_y_final - div_a_final);
    % u = u/sum(u(:));


end

