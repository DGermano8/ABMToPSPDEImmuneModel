function u = computeSimplePhenotypeModel(u, dCdx, dCdy, A, A_Ic, params)

    Dx = params.D;
    chix = params.C_chi;
    chia = params.chia;
    rho_plus = params.P_A;
    rho_minus = params.P_D;
    amax = params.activatedAge;

    dx = params.dx;
    dy = params.dy;
    da = params.da;
    dt = params.dt;

    nx = params.Nx;
    ny = params.Ny;
    na = params.Na;


    a_ind_1 = params.a_ind_1;
    a_ind_2 = params.a_ind_2;
    
    % --- Spatial Flux Term ---
    uxdiff = u(2:end, :, :) - u(1:end-1, :, :);
    uxadd = u(2:end, :, :) + u(1:end-1, :, :);

    % x-direction
    flux_west_diff = Dx * (uxdiff) / dx;
    u_west = 0.5 * (uxadd);
    flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx

    flux_west = flux_west_chem + flux_west_diff;
    flux_west = cat(1, zeros(1,ny,na), flux_west);


    flux_east_diff = Dx * (uxdiff) / dx;
    u_east = 0.5 * (uxadd); % Central difference for u at face
    flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx
    
    flux_east = flux_east_chem + flux_east_diff;
    flux_east = cat(1, flux_east, zeros(1,ny,na));

    grad_x = (flux_east - flux_west) / dx;


    % y-direction
    uydiff = u(:, 2:end, :) - u(:, 1:end-1, :);
    uyadd = u(:, 2:end, :) + u(:, 1:end-1, :);

    flux_south_diff = Dx * (uydiff) / dy;
    u_south = 0.5 * (uyadd); % Central difference for u at face
    flux_south_chem = -chix .* u_south .* dCdy; % Note the index for gradCy
    
    flux_south = flux_south_chem + flux_south_diff;
    flux_south = cat(2, zeros(nx,1,na), flux_south);


    flux_north_diff = Dx * (uydiff) / dy;
    u_north = 0.5 * (uyadd); % Central difference for u at face
    flux_north_chem = -chix .* u_north .* dCdy;   % Note the index for gradCy

    flux_north = flux_north_chem + flux_north_diff;
    flux_north  = cat(2, flux_north, zeros(nx,1,na));

    grad_y = (flux_north - flux_south) / dy;


    % --- Phenotypic Advection Term ---
    a_plus_half_a = (A(:,:,a_ind_1) + A(:,:,:)) / 2;
    f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_plus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
    f_k_plus_half_lessThanZero = (f_k_plus_half <= 0);
    u_upwind_plus = u(:, :, a_ind_1) .*(f_k_plus_half_lessThanZero) + u(:, :, :) .* (1-f_k_plus_half_lessThanZero);

    a_minus_half_a = (A(:,:,:) + A(:,:,a_ind_2) ) / 2;
    f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* A_Ic(:,:, :) - rho_minus .* a_minus_half_a./amax .* (1-A_Ic(:,:, :)) ) ;
    f_k_minus_half_lessThanZero = (f_k_minus_half <= 0);
    u_upwind_minus = u(:, :, :) .*(f_k_minus_half_lessThanZero) + u(:, :, a_ind_2) .* (1-f_k_minus_half_lessThanZero);

    flux_a_plus = chia .* f_k_plus_half .* u_upwind_plus;
    flux_a_minus = chia .* f_k_minus_half .* u_upwind_minus;

    % Boundary in a (no flux)
    flux_a_minus(:,:,1) = 0;
    % Boundary in a (no flux)
    flux_a_plus(:,:,end) = 0;

    grad_a = -( (flux_a_plus - flux_a_minus) ) / da;


    % Update u
    u = u + dt * (grad_x + grad_y + grad_a);

end