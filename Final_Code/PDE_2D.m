clear all;
% close all;

addpath(genpath('UserFunctions'))

ModelParams = struct();
rndSeed = 12;

domainBoundary = struct();
domainBoundary.x_max =  60;
domainBoundary.y_max =  60;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 100;
ModelParams.p_move = 0.25;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 16;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 3;

% T-cells
ModelParams.P_A = 0.2; % Antigen gain rate
ModelParams.P_D = 0.5;  % Antigen loss rate
ModelParams.activatedAge  = 20;

ModelParams.dx_PDE = 1/3;
ModelParams.dy_PDE = 1/3;
ModelParams.dt_PDE = 0.025;

ModelParams.plot_traj = true;
ModelParams.t_plot = 25/ModelParams.dt_PDE;


purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];

cmap = myColour3Gradient(1+ModelParams.T_final/(ModelParams.t_plot*ModelParams.dt_PDE),  yellow_colour, red_colour, purple_colour);


if ( (4*ModelParams.dt_PDE)/(ModelParams.dx_PDE.^2) > 1)

    disp('Choose smaller dx and dt for ABM')
else
    disp('Nice params!')
end
(4*ModelParams.dt_PDE)/(ModelParams.dx_PDE.^2)

%%
    
    [U, A, U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params] = PDE_SetUp(rndSeed,domainBoundary,ModelParams);
    
    if ModelParams.plot_traj
        f = figure;
        f.Position = [106 360 1360 360];
    end
   
    activatedCells = zeros(params.nt+1,1);
    tic;
    for n = 0:params.nt
        
        % disp(['time: ', num2str(n)])
        % tic;
        [U, A, U_a, A_a] = computePDEModel(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy,params);
        % toc;
        % disp(['  '])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = (A+A_a);
        % tmp = A.*U;

        activatedCells(n+1) = ( sum( tmp(:) ) );
        
        
        if mod(n, ModelParams.t_plot) == 0 && ModelParams.plot_traj
            iter = n/ModelParams.t_plot + 1;
            x = linspace(0, params.Lx, params.Nx);
            y = linspace(0, params.Ly, params.Ny);
            [X, Y] = meshgrid(x, y);

            clf;
            total_density = sum(U(:));
            title_string = ['Time = ', num2str(n * params.dt), ', total activation = ', num2str(100*activatedCells(n+1),'%.1f'), '%', ', Density = ',num2str(sum(U(:)))];
            sgtitle(title_string);
    
            subplot(1,3,1)
            surf(X, Y, (U/(sum(U(:)))),'EdgeAlpha',0);

            % hold on;
            % contour(X,Y,U, [1*10^-3, 1*10^-3], 'LineWidth',1.5, 'EdgeColor',cmap(iter,:))
            xlabel('x'); ylabel('y'); title('U');
            colorbar;
            % caxis([0 5*10^(-5)])
            view(0,90)
    
            subplot(1,3,2)
            surf(X, Y, (A+A_a)/sum(U(:)),'EdgeAlpha',0);

            % hold on;
            % contour(X,Y,A, [4*10^-4, 4*10^-4], 'LineWidth',1.5, 'EdgeColor',cmap(iter,:))

            xlabel('x'); ylabel('y'); title('A');
            colorbar; % caxis([0 1]);
            view(0,90)
    
            subplot(1,3,3)
            plot( ModelParams.dt_PDE*0:n, activatedCells(1:(n+1)), 'LineWidth',2.5)
            axis([0 params.nt 0 1])
    
            drawnow;
        end
        
    end
    toc;

% end



