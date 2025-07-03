function [U, C, dCdx, dCdy, Ic, A_Ic, Inds_Ic_st, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams)

    Lx = domainBoundary.x_max;              % Domain size in x-direction
    Ly = domainBoundary.y_max;              % Domain size in y-direction
    La = ModelParams.activatedAge;
    
    params = struct();

    dx = ModelParams.dx_PDE;
    dy = ModelParams.dy_PDE;
    da = ModelParams.da_PDE;
    
    Nx = Lx/dx + 1;               % Number of grid points in x
    Ny = Ly/dy + 1;               % Number of grid points in y
    Na = La/da + 1;
    dt = ModelParams.dt_PDE;

    params.dt = ModelParams.dt_PDE;
    
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    a = linspace(0, La, Na);

    a_ind_1 = [2:Na, Na];
    a_ind_2 = [1 1:(Na-1)];
    params.a_ind_1 = a_ind_1;
    params.a_ind_2 = a_ind_2;

    [Y,X,A] = meshgrid(y,x,a);

    nt = ceil(ModelParams.T_final / ModelParams.dt_PDE);    % Number of time steps
    params.nt = nt;
    
    U = zeros(Nx, Ny, Na);          % Population density u
    U(1,:,1) = 1;
    U=U/(sum(U(:)));

    params.dx = dx;
    params.dy = dy;
    params.da = da;
    

    params.D = ModelParams.p_move;
    params.C_chi = ModelParams.C_sens;
    params.C_scaling = ModelParams.C_scaling;
    params.P_A = ModelParams.P_A;
    params.P_D = ModelParams.P_D;
    params.Ny = Ny;
    params.Nx = Nx;
    params.Na = Na;
    params.Ly = Ly;
    params.Lx = Lx;
    params.La = La;
    params.p_move = ModelParams.p_move;
    params.activatedAge = ModelParams.activatedAge;
    activatedAge = ModelParams.activatedAge;
    
    % Initialize fields
    C = zeros(Nx, Ny, Na);          % Chemical concentration c

    R_DC = ModelParams.R_DC;

    % Get the DC clusters
    rng(rndSeed);
    if ModelParams.NumDCs > 0
        [x_DC,y_DC] = DC_Clusters(ModelParams.NumDCs,R_DC+2*R_DC,ModelParams.numberOfClusters,domainBoundary);
        for i = 1:size(x_DC, 1)
            x_DC(i) = x_DC(i) - mod(x_DC(i),dx);
            y_DC(i) = y_DC(i) - mod(y_DC(i),dy);
                
        end
        DCs = getDCField(Nx, Ny, [x_DC,y_DC], R_DC,      dx,dy); % Proximity field
        [Ic, A_Ic, Inds_Ic] =  getProximityField(Nx, Ny, Na,  [x_DC,y_DC], R_DC ,dx,dy, R_DC); % Proximity field
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

    dCdx = (C(2:end, :, :) - C(1:(end-1),:, :)) / dx;
    dCdy = (C(:, 2:end,:) - C(:,1:(end-1),:)) / dy;


    inds = [];
    inds_north = [];
    inds_south = [];
    inds_east = [];
    inds_west = [];
    for ii = 1:size(Inds_Ic)
        for ia = 1:Na
            if (Ic(Inds_Ic(ii,1)+1,Inds_Ic(ii,2)) + Ic(Inds_Ic(ii,1)-1,Inds_Ic(ii,2)) ...
                    + Ic(Inds_Ic(ii,1),Inds_Ic(ii,2)+1) + Ic(Inds_Ic(ii,1),Inds_Ic(ii,2)-1)) > 0
                inds = [inds; sub2ind([Nx, Ny, Na], Inds_Ic(ii,1),Inds_Ic(ii,2),ia)];
            end
    
            if ~(ismember([Inds_Ic(ii,1),Inds_Ic(ii,2)-1],Inds_Ic,'rows'))
                inds_north = [inds_north; sub2ind([Nx, Ny, Na], Inds_Ic(ii,1),Inds_Ic(ii,2)-1,ia)];
            end
            if ~(ismember([Inds_Ic(ii,1),Inds_Ic(ii,2)+1],Inds_Ic,'rows'))
                inds_south = [inds_south; sub2ind([Nx, Ny, Na], Inds_Ic(ii,1),Inds_Ic(ii,2)+1,ia)];
            end
            if ~(ismember([Inds_Ic(ii,1)-1,Inds_Ic(ii,2)],Inds_Ic,'rows'))
                inds_east = [inds_east; sub2ind([Nx, Ny, Na], Inds_Ic(ii,1)-1,Inds_Ic(ii,2),ia)];
            end
            if ~(ismember([Inds_Ic(ii,1)+1,Inds_Ic(ii,2)],Inds_Ic,'rows'))
                inds_west = [inds_west; sub2ind([Nx, Ny, Na], Inds_Ic(ii,1)+1,Inds_Ic(ii,2),ia)];
            end
        end
    end

    Inds_Ic_st = struct();
    Inds_Ic_st.inds = inds;
    Inds_Ic_st.inds_north = inds_north;
    Inds_Ic_st.inds_south = inds_south;
    Inds_Ic_st.inds_east = inds_east;
    Inds_Ic_st.inds_west = inds_west;


end


%%

function [Ic, A_Ic, Inds_Ic] = getProximityField(Nx, Ny, Na, DC_positions, d_threshold ,dx,dy, R_DC)
    % Compute proximity field based on distance threshold to DCs
    Ic = zeros(Nx, Ny, Na);
    A_Ic = zeros(Nx, Ny, Na);
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
                    if t_jj>Ny t_jj=Ny;end
                    if t_ii<1; t_ii=1; end
                    if t_ii>Nx t_ii=Nx;end

                    Ic(t_ii, t_jj, :) = 1;
                    Inds_Ic = [Inds_Ic; t_ii, t_jj];
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
                        if t_jj>Ny t_jj=Ny;end
                        if t_ii<1; t_ii=1; end
                        if t_ii>Nx t_ii=Nx;end
                                
                        A_Ic(t_ii, t_jj,  :) = 1;
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

