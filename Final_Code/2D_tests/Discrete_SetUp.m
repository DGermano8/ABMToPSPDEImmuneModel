function [walker_positions, walker_activation, C, DCLingerTime, DC_model, params] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams)
    
    params = struct();

    Lx = domainBoundary.x_max-domainBoundary.x_min;              % Domain size in x-direction
    Ly = domainBoundary.y_max-domainBoundary.y_min;              % Domain size in y-direction
    
    Nx = Lx/ModelParams.dx_ABM + 1;                % Number of grid points in x
    Ny = Ly/ModelParams.dy_ABM + 1;                % Number of grid points in y
    dt = ModelParams.dt_ABM;
    
    x = linspace(domainBoundary.x_min, domainBoundary.x_max, Nx);
    y = linspace(domainBoundary.y_min, domainBoundary.y_max, Ny);
    [X, Y] = meshgrid(x, y);

    nt = ceil(ModelParams.T_final / dt);    % Number of time steps
    params.nt = nt;
    
    U = zeros(Ny, Nx);          % Population density u

    dx = ModelParams.dx_ABM;
    dy = ModelParams.dy_ABM;
    da = ModelParams.da_ABM;
    params.dx = dx;
    params.dy = dy;
    params.dt = dt;
    params.da = da;

    % Initialize walker positions randomly on the lattice

    %%% For ** LEFT ** Initial Condition
    params.num_walkers = ModelParams.NumbAgents;
    walker_activation = zeros(params.num_walkers, 1);
    walker_positions = randi([1, Ny], params.num_walkers, 2); % Each row is [x, y]
    walker_positions(:,1) = 0;
    walker_positions(:,2) = dx*(walker_positions(:,2)-1);

    params.C_chi = (ModelParams.C_sens)*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2);
    params.C_scaling = ModelParams.C_scaling;

    params.P_A = (ModelParams.P_A/ModelParams.da_ABM)*ModelParams.dt_ABM;
    params.P_D = (ModelParams.P_D/ModelParams.da_ABM)*ModelParams.dt_ABM;
    params.Ny = Ny;
    params.Nx = Nx;
    params.Ly = Ly;
    params.Lx = Lx;
    params.p_move = (ModelParams.p_move)*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2);
    params.activatedAge = ModelParams.activatedAge;
    
    % Initialize fields
    C = zeros(Ny, Nx);          % Chemical concentration c
    
    R_DC = ModelParams.R_DC;

    % Get the DC clusters
    rng(rndSeed);
    if ModelParams.NumDCs > 0
        [x_DC,y_DC] = DC_Clusters(ModelParams.NumDCs,R_DC+2*R_DC,ModelParams.numberOfClusters,domainBoundary);
        for i = 1:size(x_DC, 1)
            x_DC(i) = x_DC(i) - mod(x_DC(i),dx);
            y_DC(i) = y_DC(i) - mod(y_DC(i),dy);
        end
        [Ic, DC_id] =  getProximityField(Nx, Ny,  [x_DC,y_DC], R_DC ,dx,dy , R_DC); % Proximity field
    else
        Ic = zeros(Ny, Nx);
        DC_id = zeros(Ny,Nx);
    end
    DCLingerTime = zeros(params.num_walkers, ModelParams.NumDCs);
    
    % initialise chemockine concentration
    if ModelParams.NumDCs > 0
        for ii = 1:Nx
            for jj=1:Ny
                if Ic(jj,ii) == 1

                    distance = ((X - dx*(ii-1)).^2 + (Y - dy*(jj-1)).^2);
                    C = C + exp(-2.0.*distance.*sqrt(params.C_scaling));
                end
            end
        end
    end
    
    C(Ic == 1) = 0;
    C = C/(max(C(:)));

    if ModelParams.C_sens == 0
        C = ones(Ny, Nx);
    end

    params.Cx_scale = params.dx;
    params.Cy_scale = params.dy;

    DC_model = struct();
    DC_model.Ic = Ic;
    DC_model.BoundaryDC = DC_id;
end



function [Ic, DC_id] = getProximityField(Nx, Ny, DC_positions, d_threshold ,dx,dy, R_DC)
    % Compute proximity field based on distance threshold to DCs
    Ic = zeros(Ny, Nx);
    DC_id = zeros(Ny, Nx);
    [X, Y] = meshgrid(1:Nx, 1:Ny);
    num_x = 1.0/dx;
    num_y = 1.0/dy;

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
                    iiy = i_y-floor(0.5*num_y) + jj;
                    iix = i_x-floor(0.5*num_x) + ii;
                    if iiy<1; iiy=1; end;
                    if iiy>Ny; iiy=Ny; end;
                    if iix<1; iix=1; end;
                    if iix>Nx; iix=Nx; end;
                    
                    Ic(iiy, iix) = 1;
                end
            end
        end
        
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
                        iiy = i_y-floor(0.5*num_y) + jj;
                        iix = i_x-floor(0.5*num_x) + ii;
                        if iiy<1; iiy=1; end;
                        if iiy>Ny; iiy=Ny; end;
                        if iix<1; iix=1; end;
                        if iix>Nx; iix=Nx; end;
                        
                        DC_id(iiy, iix) = i;
                    end
                end
            end
        end        
    end
end



function g = getDCField(Nx, Ny, DC_positions, d_threshold ,dx,dy)
    % Compute proximity field based on distance threshold to DCs
    g = zeros(Ny, Nx);
    [X, Y] = meshgrid(1:Nx, 1:Ny);
    for i = 1:size(DC_positions, 1)
        DC_x = DC_positions(i,1);% - mod(DC_positions(i,1),dx);
        DC_y = DC_positions(i,2);% - mod(DC_positions(i,2),dy);
        iip1 = DC_x+1; jjp1 = DC_y+1;
        iim1 = DC_x-1; jjm1 = DC_y-1;
        if iip1 > Nx; iip1=Nx; end
        if iim1 < 1; iim1=1; end
        if jjp1 > Ny; jjp1=Ny; end
        if jjm1 < 1; jjm1=1; end

        g(DC_y,DC_x) = 1;
        g(DC_y,iip1) = 1;
        g(DC_y,iim1) = 1;
        g(jjm1,DC_x) = 1;
        g(jjp1,DC_x) = 1;
        
    end
end