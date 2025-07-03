clear all;
close all;

ModelParams.dx_ABM = 1;
ModelParams.dy_ABM = 1;

domainBoundary = struct();
domainBoundary.x_max =  60;
domainBoundary.y_max =  60;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.NumDCs = 16;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 8;

Lx = domainBoundary.x_max-domainBoundary.x_min;              % Domain size in x-direction
Ly = domainBoundary.y_max-domainBoundary.y_min;

Nx = Lx/ModelParams.dx_ABM + 1;                % Number of grid points in x
Ny = Ly/ModelParams.dy_ABM + 1;
x = linspace(domainBoundary.x_min, domainBoundary.x_max, Nx);
y = linspace(domainBoundary.y_min, domainBoundary.y_max, Ny);
[X, Y] = meshgrid(x, y);
% rng(rndSeed+1);
    if ModelParams.NumDCs > 0
        [x_DC,y_DC,cluster_id] = DC_Clusters(ModelParams.NumDCs, 1+2*1,ModelParams.numberOfClusters,domainBoundary);
        for i = 1:size(x_DC, 1)
            x_DC(i) = x_DC(i) - mod(x_DC(i),ModelParams.dx_ABM);
            y_DC(i) = y_DC(i) - mod(y_DC(i),ModelParams.dx_ABM);
        end
        [Ic, DC_id] =  getProximityField(Nx, Ny,  [x_DC,y_DC], 1 ,ModelParams.dx_ABM,ModelParams.dx_ABM , 1,cluster_id); % Proximity field
    else
        Ic = zeros(Ny, Nx);
        DC_id = zeros(Ny,Nx);
    end

cmap = 1/255*[166,206,227;
31,120,180;
178,223,138;
51,160,44;
251,154,153;
227,26,28;
253,191,111;
255,127,0];

cmap = [1,1,1;
        cmap];

f = figure;
f.Position = [1513         449         298         241];
surf(X, Y,Ic,'EdgeColor','none')
colormap(cmap)
view(2)
xticks([0:20:60])
yticks([0:20:60])
axis([-0.1 60.1 -0.1 60.1])
% axis off
set(gcf,'Color',[1 1 1])

export_fig cluster_8.png -m2.5


%%
function [Ic, DC_id] = getProximityField(Nx, Ny, DC_positions, d_threshold ,dx,dy, R_DC, cluster_id)
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
                    Ic(i_y-floor(0.5*num_y) + jj, i_x-floor(0.5*num_x) + ii) = cluster_id(i);
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
                        DC_id(i_y-floor(0.5*num_y) + jj, i_x-floor(0.5*num_x) + ii) = i;
                    end
                end
            end
        end        
    end
end

