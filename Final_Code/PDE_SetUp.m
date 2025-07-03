function [U, A, U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params] = PDE_SetUp(rndSeed,domainBoundary,ModelParams)

    Lx = domainBoundary.x_max;              % Domain size in x-direction
    Ly = domainBoundary.y_max;              % Domain size in y-direction
    
    params = struct();

    dx = ModelParams.dx_PDE;
    dy = ModelParams.dy_PDE;
    
    Nx = Lx/dx + 1;               % Number of grid points in x
    Ny = Ly/dy + 1;               % Number of grid points in y
    dt = ModelParams.dt_PDE;

    params.dt = ModelParams.dt_PDE;
    
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    [X, Y] = meshgrid(x, y);

    nt = ceil(ModelParams.T_final / ModelParams.dt_PDE);    % Number of time steps
    params.nt = nt;
    
    U_a = zeros(Ny, Nx); A_a = zeros(Ny, Nx);
    U = zeros(Ny, Nx);          % Population density u
    U(:,1) = 1;
    U=U/(sum(U(:)));

    params.dx = dx;
    params.dy = dy;

    params.D = ModelParams.p_move;
    params.C_chi = ModelParams.C_sens;
    params.C_scaling = ModelParams.C_scaling;
    params.P_A = ModelParams.P_A;
    params.P_D = ModelParams.P_D;
    params.Ny = Ny;
    params.Nx = Nx;
    params.Ly = Ly;
    params.Lx = Lx;
    params.p_move = ModelParams.p_move;
    params.activatedAge = ModelParams.activatedAge;
    activatedAge = ModelParams.activatedAge;
    
    % Initialize fields
    C = zeros(Ny, Nx);          % Chemical concentration c
    A = zeros(Ny, Nx);          % Antigen level

    R_DC = ModelParams.R_DC;

    % Get the DC clusters
    rng(rndSeed);
    if ModelParams.NumDCs > 0
        [x_DC,y_DC] = DC_Clusters(ModelParams.NumDCs,R_DC+2*R_DC,ModelParams.numberOfClusters,domainBoundary);
        for i = 1:size(x_DC, 1)
            x_DC(i) = x_DC(i) - mod(x_DC(i),dx);
            y_DC(i) = y_DC(i) - mod(y_DC(i),dy);
            
            % bad_dc_x = [true true];
            % bad_dc_y = [true true];
            % while any([bad_dc_x bad_dc_y])
            %     if x_DC(i) - 1.5/dx < 0
            %         x_DC(i) = x_DC(i) + dx;
            %     else
            %         bad_dc_x(1) = false;
            %     end
            %     if x_DC(i) + 1.5/dx > Lx
            %         x_DC(i) = x_DC(i) - dx;
            %     else
            %         bad_dc_x(2) = false;
            %     end
            % 
            %     if y_DC(i) - 1.5/dy < 0
            %         y_DC(i) = y_DC(i) + dy;
            %     else
            %         bad_dc_y(1) = false;
            %     end
            %     if y_DC(i) + 1.5/dy > Ly
            %         y_DC(i) = y_DC(i) - dy;
            %     else
            %         bad_dc_y(2) = false;
            %     end
            % end

                
        end
        DCs = getDCField(Nx, Ny, [x_DC,y_DC], R_DC,      dx,dy); % Proximity field
        [Ic, A_Ic, Inds_Ic_tmp] =  getProximityField(Nx, Ny,  [x_DC,y_DC], R_DC ,dx,dy, R_DC); % Proximity field
    else
        DCs = zeros(Ny, Nx);
        Ic = zeros(Ny, Nx);
    end
    
    % initialise chemockine concentration
    if ModelParams.NumDCs > 0
        for i = 1:size(y_DC, 1)
            DC_x = x_DC(i)-0.5*dx;
            DC_y = y_DC(i)-0.5*dy;
            distance = ((X - DC_x).^2 + (Y - DC_y).^2);
            C = C + exp(-2.0.*distance.*sqrt(params.C_scaling));
        end
    else
        C = (Y-Ly)/Ly;
    end
    C(Ic == 1) = 0;
    % C = (X-Ly)/Ly;
    C = C/max(C(:));
    
    % initialise Chemotactic pde and BCs:
    C_old = C;   C_old(Ic == 1) = 0;
    C_X_old = C; C_X_old(Ic == 1) = 0;
    C_Y_old = C; C_Y_old(Ic == 1) = 0;
    
    % clf;
    % surf(Ic, 'edgecolor','none'); view(2);
    % drawnow;
    % Inds_Ic = ][;
    % for kk = 1:length(Inds_Ic_tmp)
    %     ii = Inds_Ic_tmp(kk, 1);
    %     jj = Inds_Ic_tmp(kk, 2);
    %     if Ic(ii+1,jj) == 0
    %  
    %         Inds_Ic = [Inds_Ic; ii jj];
    %     elseif Ic(ii-1,jj) == 0
    %         Inds_Ic = [Inds_Ic; ii jj];
    %     end
    % 
    %     if Ic(ii,jj+1) == 0
    %         Inds_Ic = [Inds_Ic; ii jj];
    %     elseif Ic(ii,jj-1) == 0
    %         Inds_Ic = [Inds_Ic; ii jj];
    %     end
    % end
    Inds_Ic = struct();
    Ic_vals = find(Ic == 1);
    Inds_Ic.Ic_vals = Ic_vals;
    Ic_iip1_jj = [];
    Ic_iim1_jj = [];
    Ic_ii_jjp1 = [];
    Ic_ii_jjm1 = [];

    iip1_jj = [];
    iim1_jj = [];
    ii_jjp1 = [];
    ii_jjm1 = [];
    
    sIc = size(Ic);
    % figure; imagesc(Ic)
    for kk = 1:length(Inds_Ic_tmp)
        ii = Inds_Ic_tmp(kk, 1);
        jj = Inds_Ic_tmp(kk, 2);
        if Ic(ii+1,jj) == 0
            
            Ic_iip1_jj = [Ic_iip1_jj; sub2ind(sIc,ii+1,jj)];
            iip1_jj = [iip1_jj; sub2ind(sIc,ii,jj)];

            % Inds_Ic = [Inds_Ic; ii jj];
        elseif Ic(ii-1,jj) == 0
            Ic_iim1_jj = [Ic_iim1_jj; sub2ind(sIc,ii-1,jj)];
            iim1_jj = [iim1_jj; sub2ind(sIc,ii,jj)];

            % Inds_Ic = [Inds_Ic; ii jj];
        end
        
        if Ic(ii,jj+1) == 0
            Ic_ii_jjp1 = [Ic_ii_jjp1; sub2ind(sIc,ii,jj+1)];
            ii_jjp1 = [ii_jjp1; sub2ind(sIc,ii,jj)];

            % Inds_Ic = [Inds_Ic; ii jj];
        elseif Ic(ii,jj-1) == 0
            Ic_ii_jjm1 = [Ic_ii_jjm1; sub2ind(sIc,ii,jj-1)];
            ii_jjm1 = [ii_jjm1; sub2ind(sIc,ii,jj)];

            % Inds_Ic = [Inds_Ic; ii jj];
        end
    end
    Inds_Ic.Ic_iip1_jj = Ic_iip1_jj;
    Inds_Ic.iip1_jj = iip1_jj;

    Inds_Ic.Ic_iim1_jj = Ic_iim1_jj;
    Inds_Ic.iim1_jj = iim1_jj;

    Inds_Ic.Ic_ii_jjp1 = Ic_ii_jjp1;
    Inds_Ic.ii_jjp1 = ii_jjp1;

    Inds_Ic.Ic_ii_jjm1 = Ic_ii_jjm1;
    Inds_Ic.ii_jjm1 = ii_jjm1;

    for ii=2:Ny-1
        for jj=2:Nx-1
            Iter_X = [ii+1; ii-1]; Iter_Y = [jj-1; jj+1];
            if Ic(ii,jj) == 1
                if Ic(ii+1,jj) == 0
                    C_X_old(ii,jj) =  C(ii+1,jj);
                    if Ic(ii-1,jj) == 1
                        C_X_old(ii-1,jj) =  C(ii+1,jj);
                    end
                end
                if Ic(ii-1,jj) == 0
                    C_X_old(ii,jj) =  C(ii-1,jj);
                    if Ic(ii+1,jj) == 1
                        C_X_old(ii+1,jj) =  C(ii-1,jj);
                    end
                end
                if Ic(ii,jj+1) == 0
                    C_Y_old(ii,jj) =  C(ii,jj+1);
                    if Ic(ii,jj-1) == 1
                        C_Y_old(ii,jj-1) =  C(ii,jj+1);
                    end
                end
                if Ic(ii,jj-1) == 0
                    C_Y_old(ii,jj) =  C(ii,jj-1);
                    if Ic(ii,jj+1) == 1
                        C_Y_old(ii,jj+1) =  C(ii,jj-1);
                    end
                end
            end
        end
    end
    c_X_extended = padarray(C_X_old, [1, 1], 'replicate'); % Extend boundaries by mirroring
    c_Y_extended = padarray(C_Y_old, [1, 1], 'replicate'); % Extend boundaries by mirroring
    % Extract shifted matrices (inner region only)
    c_shift_xp = c_X_extended(3:end, 2:end-1); % Shift down
    c_shift_xm = c_X_extended(1:end-2, 2:end-1); % Shift up
    c_shift_yp = c_Y_extended(2:end-1, 3:end); % Shift right
    c_shift_ym = c_Y_extended(2:end-1, 1:end-2); % Shift left
    % For chemotactic gradients
    dc1_dx = (c_shift_xp - C_X_old) / params.dx;
    dc2_dx = (C_X_old - c_shift_xm) / params.dx;
    dc1_dy = (c_shift_yp - C_Y_old) / params.dy;
    dc2_dy = (C_Y_old - c_shift_ym) / params.dy;

    U(Ic == 1) = 0;
    A(Ic == 1) = 0;
    A_a(Ic == 1) = 0;

    % find(Ic == 1);
end


%%

function [Ic, A_Ic, Inds_Ic] = getProximityField(Nx, Ny, DC_positions, d_threshold ,dx,dy, R_DC)
    % Compute proximity field based on distance threshold to DCs
    Ic = zeros(Ny, Nx);
    A_Ic = zeros(Ny, Nx);
    [X, Y] = meshgrid(1:Nx, 1:Ny);
    num_x = 1.0/dx;
    num_y = 1.0/dy;
    Inds_Ic = [];

    for i = 1:size(DC_positions, 1)
        DC_x = DC_positions(i,1);
        DC_y = DC_positions(i,2);

        x_id = round(DC_x/dx);
        y_id = round(DC_y/dy);

        Iter_X = [x_id; x_id+num_x; x_id-num_x; x_id; x_id];
        Iter_Y = [y_id; y_id; y_id; y_id+num_y; y_id-num_y];
        % center
        for kk=1:length(Iter_X)
            i_x = Iter_X(kk);
            i_y = Iter_Y(kk);
            for ii=1:num_x
                for jj=1:num_y
                    t_jj = i_y-floor(0.5*num_y) + jj;
                    t_ii = i_x-floor(0.5*num_x) + ii;
                    if t_jj<1; t_jj=1; end
                    if t_jj>Nx t_jj=Nx;end
                    if t_ii<1; t_ii=1; end
                    if t_ii>Ny t_ii=Ny;end

                    Ic(t_jj, t_ii) = 1;
                    Inds_Ic = [Inds_Ic; t_jj, t_ii];
                end
            end
        end
        % disp('done')
        for kk=1:length(Iter_X)
            xx_id = Iter_X(kk);
            yy_id = Iter_Y(kk);
            
            Iter_XX = [xx_id; xx_id+num_x; xx_id-num_x; xx_id;       xx_id;       xx_id-num_x; xx_id-num_x; xx_id+num_x; xx_id+num_x];
            Iter_YY = [yy_id; yy_id;       yy_id;       yy_id+num_y; yy_id-num_y; yy_id-num_y; yy_id+num_y; yy_id-num_y; yy_id+num_y];
            % center
            for kk=1:length(Iter_XX)
                i_x = Iter_XX(kk);
                i_y = Iter_YY(kk);
                for ii=1:num_x
                    for jj=1:num_y
                        t_jj = i_y-floor(0.5*num_y) + jj;
                        t_ii = i_x-floor(0.5*num_x) + ii;
                        if t_jj<1; t_jj=1; end
                        if t_jj>Nx t_jj=Nx;end
                        if t_ii<1; t_ii=1; end
                        if t_ii>Ny t_ii=Ny;end
                                
                        A_Ic(t_jj, t_ii) = 1;
                    end
                end
            end
        end        
        % g(distance < d_threshold) = 1;
    end
    % A_Ic = A_Ic - Ic;
end


function g = getDCField(Nx, Ny, DC_positions, d_threshold ,dx,dy)
    % Compute proximity field based on distance threshold to DCs
    g = zeros(Ny, Nx);
    [X, Y] = meshgrid(1:Nx, 1:Ny);
    for i = 1:size(DC_positions, 1)
        DC_x = DC_positions(i,1) - mod(DC_positions(i,1),dx);
        DC_y = DC_positions(i,2) - mod(DC_positions(i,2),dy);
        distance = sqrt((X*dx - DC_x).^2 + (Y*dy - DC_y).^2);
        g(distance <= d_threshold) = 1;
    end
end

