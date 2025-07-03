clc;
% close all;
clear all;

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];


%%

ModelParams = struct();

domainBoundary = struct();
domainBoundary.x_max =  50;
domainBoundary.y_max =  50;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 50;
% ModelParams.T_final = 1;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 1;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 0.5; % Antigen gain rate
ModelParams.P_D = 0.5;  % Antigen loss rate
ModelParams.activatedAge  = 50;

ModelParams.dx_ABM = 0.5;
ModelParams.dy_ABM = 0.5;
ModelParams.dt_ABM = 0.1;            % Time step
ModelParams.da_ABM = 1;            % Time step
ModelParams.NumbAgents = 10000;

ModelParams.dx_PDE = 0.5;
ModelParams.dy_PDE = 0.5;
ModelParams.dt_PDE = 0.1;
ModelParams.da_PDE = 1;

ModelParams.plot_traj = true;
recordRate = 800;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;


rndSeed = 4;


poolobj = gcp("nocreate"); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(6);
end

P_A = [0:0.25:1];
P_D = [0:0.25:1];
chia_x = [0.5];
% chia_x = [0.5];
Diff_x = [0.5];

ABM_Activation = zeros(length(chia_x),length(Diff_x),length(P_A),length(P_D));
PDE_Activation = zeros(length(chia_x),length(Diff_x),length(P_A),length(P_D));

iter = 0;
tic;


for Xx = 1:length(chia_x)
for Dx = 1:length(Diff_x)
    ModelParams.C_sens = chia_x(Xx);
    ModelParams.p_move = Diff_x(Dx);

    for pA = 1:length(P_A);
    
        ModelParams.P_A = P_A(pA); % Antigen gain rate  
    
        parfor pD = 1:length(P_D);
            % [pA, pD]
            
            % ModelParams.P_D = P_D(pD);  % Antigen loss rate
        
            [u, ~, dCdx, dCdy, Ic, A_Ic, Inds_Ic_st, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams);
            activatedCells_PDE = zeros(params.nt,1);
        
            [walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);    
            activatedCells_ABM = zeros(params_ABM.nt,1);
    
            params.P_D = P_D(pD);
            params_ABM.P_D = P_D(pD);
    
            x = linspace(0, params.Lx, params.Nx);
            y = linspace(0, params.Ly, params.Ny);
            a = linspace(0, params.La, params.Na);
        
            for n=0:params.nt
                u = computePhenotypeModel(u, dCdx, dCdy, A, Inds_Ic_st, A_Ic, params);
                activatedCells_PDE(n+1) = squeeze(sum(u,[1,2]))'*a'/params.activatedAge;
                
                [walker_positions, walker_activation, C, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
                activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*params_ABM.num_walkers);
            end
            
            ABM_Activation(Xx,Dx,pA,pD) = activatedCells_ABM(end);
            PDE_Activation(Xx,Dx,pA,pD) = activatedCells_PDE(end);
            
        end
    
    
        iter = iter + length(P_D);
        prog = (iter)/(length(chia_x)*length(Diff_x)*length(P_D)*length(P_A));
        elapsedTime = toc;
        estimatedTotalTime = elapsedTime / prog;
        remainingTime = (estimatedTotalTime - elapsedTime)/60;
    
        % disp(sprintf('\rProgress: %1.2f%%, Estimated time remaining: %1.2f minutes', 100*prog,  remainingTime))
        clc;
        fprintf('Progress: %1.2f%%, Estimated time remaining: %1.2f minutes', 100*prog,  remainingTime);
        fprintf('\n');
    
    end
end
end
disp('\n');
toc;


%%


[pA, pD] = meshgrid(P_A,P_D);

C_levels = [0:0.1:1];

cmap = myColour3Gradient(length(C_levels)+1, yellow_colour, red_colour, purple_colour);

figure;

iter = 0;
for Dx = 1:length(Diff_x)
for Xx = 1:length(chia_x)
    iter = iter + 1;
    ABM_XD = reshape(ABM_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
    PDE_XD = reshape(PDE_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);

    subplot(length(Diff_x),length(chia_x),iter)
    hold on;

    contour(pA, pD, ABM_XD', C_levels, ...
        'linestyle', '--',  'LineWidth', 3, ...
        'ShowText','on','LabelColor','w')
    colormap(cmap);

    contour(pA, pD, PDE_XD', C_levels, ...
        'linestyle', '--', 'LineWidth', 3 , 'color', back_colour, ...
        'ShowText','on','LabelColor','k')

    % surf(pA, pD, ABM_XD')
    % colormap(cmap);
    % surf(pA, pD, PDE_XD'-ABM_XD')
    % colormap(cmap);

    hold off;

    % title(['$D = ', num2str(Diff_x(Dx)), ', \chi = ' num2str(chia_x(Xx)) '$'] , 'FontSize',14 , 'Interpreter','latex')
    title(['$\chi = ' num2str(chia_x(Xx)) '$'] , 'FontSize',14 , 'Interpreter','latex')
    
    xlabel('$\mu_+$', 'FontSize',16, 'Interpreter','latex')
    ylabel('$\mu_-$', 'FontSize',16, 'Interpreter','latex')
    c = colorbar;

    c.Ticks = C_levels;
    c.Label.String = 'Proportion of T cell activation';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 12;

    caxis([0 1])

end
end
%%
f = gcf;
f.Color = [1 1 1];
% f.Position = [1728 381 1215 250];
% % export_fig 2D_SingleDC.png -m2.5
% 
    

%%






[pA, pD] = meshgrid(P_A,P_D);
C_levels = [0:0.1:1];
cmap = myColour3Gradient(255, yellow_colour, red_colour, purple_colour);
figure;

iter = 0;
for Dx = 1:length(Diff_x)
for Xx = 1:length(chia_x)
    iter = iter + 1;
    ABM_XD = reshape(ABM_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
    PDE_XD = reshape(PDE_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);

    subplot(1,3,1)
    surf(pA, pD, ABM_XD', 'EdgeColor','none')
    colormap(cmap);
    view(2);
    c = colorbar;
    c.Ticks = C_levels;
    caxis([0 1])

    subplot(1,3,2)
    surf(pA, pD, PDE_XD', 'EdgeColor','none')
    view(2);
    c = colorbar;
    c.Ticks = C_levels;
    caxis([0 1])

    subplot(1,3,3)
    surf(pA, pD, ABM_XD'-PDE_XD', 'EdgeColor','none')
    colormap(cmap);
    view(2);
    c = colorbar;
    
end
end
















