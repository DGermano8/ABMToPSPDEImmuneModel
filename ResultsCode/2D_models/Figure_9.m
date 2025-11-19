
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
domainBoundary.x_max =  10;
domainBoundary.y_max =  10;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 500;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 1.0;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 1;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 0.25; % Antigen gain rate
ModelParams.P_D = 0.9;  % Antigen loss rate
ModelParams.activatedAge  = 50;

dx = 0.5;
dy = dx;
dt = 0.05;
da = 0.1;

ModelParams.dx_ABM = dx;
ModelParams.dy_ABM = dy;
ModelParams.dt_ABM = dt;            % Time step
ModelParams.da_ABM = da;            % Time step
ModelParams.NumbAgents = 5000;

ModelParams.dx_PDE = dx;
ModelParams.dy_PDE = dy;
ModelParams.dt_PDE = dt;
ModelParams.da_PDE = da;

ModelParams.plot_traj = true;
recordRate = 20/dt;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;

rndSeed = 1;

%%


poolobj = gcp("nocreate"); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(6);
end

P_A = [0.01 0.05 0.1:0.1:1];
P_D = [0.01 0.05 0.1:0.1:1];


chia_x = [0.5];
% chia_x = [0.5];
Diff_x = [0.5];
a_max = [20 50 100];

ABM_Activation = zeros(length(a_max),length(Diff_x),length(P_A),length(P_D));
PDE_Activation = zeros(length(a_max),length(Diff_x),length(P_A),length(P_D));
APP_Activation = zeros(length(a_max),length(Diff_x),length(P_A),length(P_D));


iter = 0;
tic;
for Xx = 1:length(a_max)
    ModelParams.da_ABM = a_max(Xx)/200;
    ModelParams.da_PDE = a_max(Xx)/200;
    for Dx = 1:length(Diff_x)


        ModelParams.C_sens = chia_x(1);
        ModelParams.p_move = Diff_x(Dx);

        ModelParams.activatedAge = a_max(Xx);
    
        for pA = 1:length(P_A);
        
            ModelParams.P_A = P_A(pA); % Antigen gain rate  
        
            parfor pD = 1:length(P_D);

                [ABM_i, PDE_i, APP_i] = build_ABM_PDE_models_and_simulate(rndSeed,domainBoundary,ModelParams, P_D(pD) );
                ABM_Activation(Xx,Dx,pA,pD) = ABM_i;
                PDE_Activation(Xx,Dx,pA,pD) = PDE_i;
                APP_Activation(Xx,Dx,pA,pD) = APP_i;

            end

            iter = iter + length(P_D);
            prog = (iter)/(length(a_max)*length(Diff_x)*length(P_D)*length(P_A));
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


%%



[pA, pD] = meshgrid(P_A,P_D);

C_levels = [0:0.1:1];

cmap = myColour3Gradient(length(C_levels)+1, yellow_colour, red_colour, purple_colour);

figure;

iter = 0;
for Dx = 1:length(Diff_x)
    for Xx = 1:length(a_max)
        iter = iter + 1;
        ABM_XD = reshape(ABM_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        PDE_XD = reshape(PDE_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        APP_XD = reshape(APP_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);

        
    
        subplot(length(Diff_x),length(a_max),iter)
        hold on;

        contour(pA, pD, APP_XD', C_levels, ...
            'linestyle', ':', 'LineWidth', 3 , 'color', blue_colour, ...
            'ShowText','on','LabelColor','w')
    
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
        title(['$\chi = ' num2str(a_max(Xx)) '$', ', $D = ' num2str(Diff_x(Dx)) '$'] , 'FontSize',14 , 'Interpreter','latex')
        
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

figure;

load('results/1DC_Sweep/Tf_100.mat')
[pA, pD] = meshgrid(P_A,P_D);
C_levels = [0:0.1:1];
cmap = myColour3Gradient(length(C_levels)+1, yellow_colour, red_colour, purple_colour);

iter_i = 1;
iter = 0;
for Dx = 1:length(Diff_x)
    for Xx = 1:length(a_max)
        ABM_XD = reshape(ABM_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        PDE_XD = reshape(PDE_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        APP_XD = reshape(APP_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);

        subplot(length(a_max),4,iter_i+4*iter)
        hold on;
        contour(pA, pD, APP_XD', C_levels, 'linestyle', ':', 'LineWidth', 2.5 , 'color', blue_colour, 'ShowText','off', 'EdgeAlpha',0.4)

        contour(pA, pD, ABM_XD', C_levels,'linestyle', '--',  'LineWidth', 3, 'ShowText','on','LabelColor','w')
        colormap(cmap);
   
        contour(pA, pD, PDE_XD', C_levels, 'linestyle', '--', 'LineWidth', 3 , 'color', back_colour, 'ShowText','on','LabelColor','k')
        hold off;
        axis tight
    
        ylabel('$\mu_-$', 'FontSize',16, 'Interpreter','latex')

        text(gca, -0.4, 0.5, ['$A_{max}=', num2str(a_max(Xx)),'$'], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
             'FontSize',16, 'Interpreter','latex');

        if iter==2
            xlabel('$\mu_+$', 'FontSize',16, 'Interpreter','latex')
        end
        if iter==0
            title('Time = 100' , 'FontSize',14 , 'Interpreter','latex')
        end

        iter = iter + 1;
    end
end

clear all
load('results/1DC_Sweep/Tf_500.mat')

iter_i = 2;
iter = 0;

for Dx = 1:length(Diff_x)
    for Xx = 1:length(a_max)
        ABM_XD = reshape(ABM_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        PDE_XD = reshape(PDE_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        APP_XD = reshape(APP_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);

        subplot(length(a_max),4,iter_i+4*iter)
        hold on;
        contour(pA, pD, APP_XD', C_levels, 'linestyle', ':', 'LineWidth', 2.5 , 'color', blue_colour, 'ShowText','off', 'EdgeAlpha',0.4)

        contour(pA, pD, ABM_XD', C_levels,'linestyle', '--',  'LineWidth', 3, 'ShowText','on','LabelColor','w')
        colormap(cmap);
   
        contour(pA, pD, PDE_XD', C_levels, 'linestyle', '--', 'LineWidth', 3 , 'color', back_colour, 'ShowText','on','LabelColor','k')
        hold off;
        axis tight
    
       
        if iter==2
            xlabel('$\mu_+$', 'FontSize',16, 'Interpreter','latex')
        end
        if iter==0
            title('Time = 500' , 'FontSize',14 , 'Interpreter','latex')
        end
        iter = iter + 1;
    end
end


clear all
load('results/1DC_Sweep/Tf_1000.mat')

iter_i = 3;
iter = 0;

for Dx = 1:length(Diff_x)
    for Xx = 1:length(a_max)
        ABM_XD = reshape(ABM_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        PDE_XD = reshape(PDE_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        APP_XD = reshape(APP_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);

        subplot(length(a_max),4,iter_i+4*iter)
        hold on;
        contour(pA, pD, APP_XD', C_levels, 'linestyle', ':', 'LineWidth', 2.5 , 'color', blue_colour, 'ShowText','off', 'EdgeAlpha',0.4)

        contour(pA, pD, ABM_XD', C_levels,'linestyle', '--',  'LineWidth', 3, 'ShowText','on','LabelColor','w')
        colormap(cmap);
   
        contour(pA, pD, PDE_XD', C_levels, 'linestyle', '--', 'LineWidth', 3 , 'color', back_colour, 'ShowText','on','LabelColor','k')
        hold off;
        axis tight


        if iter==2
            xlabel('$\mu_+$', 'FontSize',16, 'Interpreter','latex')
        end
        if iter==0
            title('Time = 1000' , 'FontSize',14 , 'Interpreter','latex')
        end
        
        iter = iter + 1;
    end
end


clear all
load('results/1DC_Sweep/Tf_5000.mat')

iter_i = 4;
iter = 0;

for Dx = 1:length(Diff_x)
    for Xx = 1:length(a_max)
        ABM_XD = reshape(ABM_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        PDE_XD = reshape(PDE_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);
        APP_XD = reshape(APP_Activation(Xx,Dx,:,:), [length(P_A),length(P_D)]);

        subplot(length(a_max),4,iter_i+4*iter)
        hold on;
        contour(pA, pD, APP_XD', C_levels, 'linestyle', ':', 'LineWidth', 3 , 'color', blue_colour, 'ShowText','on','LabelColor','w' , 'EdgeAlpha',1)

        contour(pA, pD, ABM_XD', C_levels,'linestyle', '--',  'LineWidth', 3, 'ShowText','on','LabelColor','w')
        colormap(cmap);
   
        contour(pA, pD, PDE_XD', C_levels, 'linestyle', '--', 'LineWidth', 3 , 'color', back_colour, 'ShowText','on','LabelColor','k')
        hold off;
        axis tight
    
        % title(['$\chi = ' num2str(a_max(Xx)) '$', ', $D = ' num2str(Diff_x(Dx)) '$'] , 'FontSize',14 , 'Interpreter','latex')
        
        if iter==2
            xlabel('$\mu_+$', 'FontSize',16, 'Interpreter','latex')
        end
        if iter==0
            title('Time = 5000' , 'FontSize',14 , 'Interpreter','latex')
        end
        % ylabel('$\mu_-$', 'FontSize',16, 'Interpreter','latex')
        c = colorbar;
    
        c.Ticks = C_levels(1:2:end);
        c.Label.String = 'T cell activation';
        c.Label.Interpreter = 'latex';
        c.Label.FontSize = 14;
        caxis([0 1])
        % c.Layout.Tile = 'east';
        
        
        iter = iter + 1;
    end
end

ax1 = subplot(3,4,1);

for ii=[4 8 12]
    ax = subplot(3,4,ii);
    ax.Position(3) = ax1.Position(3);
end


%%

f = gcf;
f.Color = [1 1 1];
% export_fig 1DC_Actication_Time.png -m2.5

