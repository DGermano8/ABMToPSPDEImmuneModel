
clear all;
%%

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];
font_colour = 1/255*[217 217 217];

%%

ModelParams = struct();

domainBoundary = struct();
domainBoundary.x_max =  150;
domainBoundary.y_max =  150;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

NumberOfTimeSamples = 2*2*7;
ModelParams.T_final = 2*24*60; 7*24*60;
ModelParams.p_move = 50 * 1/16;  % t is measured in seconds

% ChemoTaxis params
ModelParams.C_sens = 1.0*ModelParams.p_move;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 128;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 1; % Antigen gain rate
ModelParams.P_D = 1;  % Antigen loss rate
ModelParams.activatedAge  = 100;

% dx = 0.2;  
% dy = dx;
% dt = 0.01;
% da = 0.01;
dx = 0.5;
dy = dx;
dt = 0.01;
da = 0.25;

ModelParams.dx_ABM = dx;
ModelParams.dy_ABM = dy;
ModelParams.dt_ABM = dt;            % Time step
ModelParams.da_ABM = da;            % Time step
ModelParams.NumbAgents = 1000;

ModelParams.dx_PDE = dx;
ModelParams.dy_PDE = dy;
ModelParams.dt_PDE = dt;
ModelParams.da_PDE = da;

ModelParams.plot_traj = true;
recordRate = 20/dt;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;

[(ModelParams.p_move)*(4*dt)/(dx.^2) (ModelParams.C_sens)*(4*dt)/(dx.^2) (ModelParams.P_A/da)*dt (ModelParams.P_D/da)*dt]
%%


poolobj = gcp("nocreate"); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(5);
end

tst_ps = [0.01 0.1 1];
P_A = tst_ps;
P_D = tst_ps;
NumberOfClusters = 2.^[0:log2(ModelParams.NumDCs)];
NumbRepeats = 5;

ABM_NumbDCAct = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,ModelParams.NumDCs+1, NumberOfTimeSamples);


ABM_Activation = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);
PDE_Activation = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);
APP_Activation = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);

ABM_FullActivation = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);
PDE_FullActivation = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);
APP_FullActivation = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);

ABM_Heterogeneity_mean = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);
ABM_Heterogeneity_var = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);
ABM_TCellAroundDC_mean = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);
ABM_TCellAroundDC_var = zeros(length(P_A),length(P_D),length(NumberOfClusters),NumbRepeats,NumberOfTimeSamples);


iter = 0;
tic;


for pA = 1:length(P_A);
for pD = 1:length(P_D);
for clst = 1:length(NumberOfClusters);

    parfor rndSeed = 1:NumbRepeats;

        ModelParams_i = struct();
        ModelParams_i = ModelParams;
        
        ModelParams_i.numberOfClusters = NumberOfClusters(clst);
        ModelParams_i.P_A = P_A(pA); % Antigen gain rate  
        ModelParams_i.P_D = P_D(pD); % Antigen gain rate  

        [ABM_i, PDE_i, APP_i, ABM_DCdist] = build_ABM_PDE_models_and_simulate_clustering_2(rndSeed,domainBoundary,ModelParams_i,NumberOfTimeSamples);
        ABM_NumbDCAct(pA,pD,clst,rndSeed,:,:) = ABM_DCdist(:,:);

        ABM_Activation(pA,pD,clst,rndSeed,:) = ABM_i(1,:);
        PDE_Activation(pA,pD,clst,rndSeed,:) = PDE_i(1,:);
        APP_Activation(pA,pD,clst,rndSeed,:) = APP_i(1,:);

        ABM_FullActivation(pA,pD,clst,rndSeed,:) = ABM_i(2,:);
        PDE_FullActivation(pA,pD,clst,rndSeed,:) = PDE_i(2,:);
        APP_FullActivation(pA,pD,clst,rndSeed,:) = APP_i(2,:);

        ABM_Heterogeneity_mean(pA,pD,clst,rndSeed,:) = ABM_i(3,:);
        ABM_Heterogeneity_var(pA,pD,clst,rndSeed,:) = ABM_i(4,:);

        ABM_TCellAroundDC_mean(pA,pD,clst,rndSeed,:) = ABM_i(5,:);
        ABM_TCellAroundDC_var(pA,pD,clst,rndSeed,:) = ABM_i(6,:);
    end

    iter = iter + NumbRepeats;
    prog = (iter)/(length(NumberOfClusters)*length(P_A)*length(P_D)*NumbRepeats);
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

simulated_time = toc

%%

figure;

times_to_plot = [NumberOfTimeSamples];
marker_style = ["o-", "square:", "^:"];
cmap = myColour3Gradient(3, yellow_colour, red_colour, purple_colour);


for its = 1:length(times_to_plot)
    
   
    av_ABM_Activation = mean(ABM_Activation(:,:,:,:,times_to_plot(its)),4);
    av_PDE_Activation = mean(PDE_Activation(:,:,:,:,times_to_plot(its)),4);
    av_APP_Activation = mean(APP_Activation(:,:,:,:,times_to_plot(its)),4);
    
    av_ABM_FullActivation = mean(ABM_FullActivation(:,:,:,:,times_to_plot(its)),4);
    av_PDE_FullActivation = mean(PDE_FullActivation(:,:,:,:,times_to_plot(its)),4);
    av_APP_FullActivation = mean(APP_FullActivation(:,:,:,:,times_to_plot(its)),4);

    mean_ABM_Activation = mean(ABM_Heterogeneity_mean(:,:,:,:,times_to_plot(its)),4);
    var_ABM_Activation = mean(ABM_Heterogeneity_var(:,:,:,:,times_to_plot(its)),4);
    

    mean_TCellAroundDC = mean(ABM_TCellAroundDC_mean(:,:,:,:,times_to_plot(its)),4);
    var_TCellAroundDC = mean(ABM_TCellAroundDC_var(:,:,:,:,times_to_plot(its)),4);
    
    
    C_levels = [0:0.1:1];
        
    
    Cells_Pre_Cluster = ModelParams.NumDCs./NumberOfClusters;
    
    
    for pA = 1:length(P_A)
    for pD = 1:length(P_D)


        marker_i = marker_style(1);
        
        color_i = pD;
        subplot(2,length(P_A),pA)
        
        hold on;
        plot(log2(Cells_Pre_Cluster) , squeeze(av_ABM_Activation(pA,pD,:)), marker_i, "Color", cmap(color_i,:), 'LineWidth',1)
        % plot(log2(NumberOfClusters) , squeeze(av_PDE_Activation(pA,pD,:)), 'o-', "Color", yellow_colour, 'LineWidth',2)
        plot(log2(Cells_Pre_Cluster) , squeeze(av_APP_Activation(pA,pD,:)), '--', "Color", cmap(color_i,:), 'LineWidth',2)
        hold off;
        axis([min(log2(Cells_Pre_Cluster))-0.5 max(log2(Cells_Pre_Cluster))+0.5 0 1])
        title(['$\mu_+$ = ', num2str(P_A(pA))], 'Interpreter','latex', 'FontSize',14)
    
        xlabel('$\log_2$(Cluster Size)', 'Interpreter','latex', 'FontSize',14)
    
    end
    end
    subplot(2,length(P_A),1)
    hold on;
    ylabel('Prop. T cell activation', 'Interpreter','latex', 'FontSize',14)
    
    % figure;
    
    Cells_Pre_Cluster = ModelParams.NumDCs./NumberOfClusters;
    
    iter = 0;
    for pA = 1:length(P_A)
    for pD = 1:length(P_D)
        
        av_DC_dist = squeeze(mean(ABM_NumbDCAct(pA,pD,:,:,:,times_to_plot(its)),4));

        % mean_DC_ints = zeros(length(Cells_Pre_Cluster),1);
        % 
        % for jj=1:length(Cells_Pre_Cluster)
        %     % min_DC = min(find(av_DC_dist(jj,:) ~= 0));
        %     % if min_DC > 1
        %     %     min_DC = max_DC - min_DC - 1;
        %     % end
        %     % max_DC = max(find(av_DC_dist(jj,:) ~= 0));
        %     % if max_DC < ModelParams.NumDCs
        %     %     max_DC = max_DC + 1;
        %     % end
        %     % 
        %     % x_fil = log2(Cells_Pre_Cluster(jj))*ones(1,(max_DC - min_DC)+1);
        %     % y_fil = min_DC:max_DC;
        %     % z_fil = [av_DC_dist(jj,min_DC:max_DC)];
        % 
        %     x_fil = log2(Cells_Pre_Cluster(jj))*ones(1,ModelParams.NumDCs+1);
        %     y_fil = 0:ModelParams.NumDCs;
        %     z_fil = [av_DC_dist(jj,:)];
        % 
        %     mean_DC_ints(jj) = z_fil*y_fil';
        % 
        % end


        marker_i = marker_style(1);

        color_i = pD;
        subplot(2,length(P_A),pA+length(P_A))
        hold on;
        plot(log2(Cells_Pre_Cluster) , squeeze(mean_ABM_Activation(pA,pD,:)), marker_i, "Color", cmap(color_i,:), 'LineWidth',1)
        % plot(log2(Cells_Pre_Cluster) , mean_DC_ints, '*', "Color", cmap(color_i,:), 'LineWidth',1)
        
        plot(log2(Cells_Pre_Cluster) , squeeze(mean_ABM_Activation(pA,pD,:)) - squeeze(var_ABM_Activation(pA,pD,:)), '--', "Color", cmap(color_i,:), 'LineWidth',1)
        plot(log2(Cells_Pre_Cluster) , squeeze(mean_ABM_Activation(pA,pD,:)) + squeeze(var_ABM_Activation(pA,pD,:)), '--', "Color", cmap(color_i,:), 'LineWidth',1)
        
        % cmap2 = myColour3Gradient(length(Cells_Pre_Cluster), yellow_colour, red_colour, purple_colour);
        % for jj=1:length(Cells_Pre_Cluster)
        %     x_fil = log2(Cells_Pre_Cluster(jj))*ones(1,ModelParams.NumDCs);
        %     y_fil = 1:ModelParams.NumDCs;
        %     z_fil = av_DC_dist(jj,:);
        %     fill3([x_fil fliplr(x_fil)], [y_fil zeros(size(y_fil))], [z_fil z_fil], ...
        %     cmap2(jj,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
        % 
        %     plot3(x_fil, y_fil, z_fil, '-', 'LineWidth',2, 'Color', cmap2(jj,:))
        % end

       
        hold off;
        axis([min(log2(Cells_Pre_Cluster))-0.5 max(log2(Cells_Pre_Cluster))+0.5 0 ModelParams.NumDCs])
        % title(['pA = ', num2str(P_A(pA))], 'Interpreter','latex', 'FontSize',14)
        title(['$\mu_+$ = ', num2str(P_A(pA))], 'Interpreter','latex', 'FontSize',14)

        xlabel('$\log_2$(Cluster Size)', 'Interpreter','latex', 'FontSize',14)
    
    end
    end
    subplot(2,length(P_A),length(P_A)+1)
    hold on;
    ylabel('\# DCs to activate', 'Interpreter','latex', 'FontSize',14)
    

    
    % figure;
    
    % Cells_Pre_Cluster = ModelParams.NumDCs./NumberOfClusters;
    % 
    % iter = 0;
    % for pA = 1:length(P_A)
    % for pD = 1:length(P_D)
    % 
    %     color_i = pD;
    %     subplot(length(P_A),length(P_D),pA+2*length(P_A))
    %     hold on;
    %     plot(log2(Cells_Pre_Cluster) , squeeze(mean_TCellAroundDC(pA,pD,:)), 'o-', "Color", cmap(color_i,:), 'LineWidth',2)
    %     % plot(log2(Cells_Pre_Cluster) , squeeze(mean_TCellAroundDC(pA,pD,:)) - squeeze(var_TCellAroundDC(pA,pD,:)), '--', "Color", cmap(color_i,:), 'LineWidth',1)
    %     % plot(log2(Cells_Pre_Cluster) , squeeze(mean_TCellAroundDC(pA,pD,:)) + squeeze(var_TCellAroundDC(pA,pD,:)), '--', "Color", cmap(color_i,:), 'LineWidth',1)
    % 
    %     hold off;
    %     axis([min(log2(Cells_Pre_Cluster))-1 max(log2(Cells_Pre_Cluster))+1 0 ModelParams.T_final])
    %     title(['pA = ', num2str(P_A(pA))])
    % 
    % end
    % end
    % subplot(length(P_A),length(P_D),2*length(P_A)+1)
    % hold on;
    % ylabel('Number T cells Per DC')
    
    
    formatSpec = '%.2f';
    
    subplot(2,length(P_A),1)
    hold on;
    p = [];
    str = [];
    for pD = 1:length(P_D)
        p = [p; plot(-1 , -1, 'o-', "Color", cmap(pD,:), 'LineWidth',2)];
        str = [str; "$\mu_-$ = " + num2str(P_D(pD),formatSpec)];
    end
    legend(p, str,'Location','northwest', 'Interpreter','latex', 'FontSize',12)

end

f = gcf;
f.Color = [1 1 1];
f.Position = [252 216 890 466];
%%
ax = subplot(2,3,6)
exportgraphics(ax,'Figs/MainResultsFig_Full_subplot_6.png','Resolution',500)

%%

f = gcf;
f.Color = [1 1 1];
f.Position = [252 216 890 466];
% export_fig MainResultsFig_Full.png -m2.5
% MakeDark()



%%
%
%
%%

figure;

times_to_plot = [4:4:NumberOfTimeSamples];
marker_style = ["o", "square:", "^:"];
cmap = myColour3Gradient(length(times_to_plot), yellow_colour, red_colour, purple_colour);


for its = 1:length(times_to_plot)
    
    
   
    av_ABM_Activation = mean(ABM_Activation(:,:,:,:,times_to_plot(its)),4);
    av_PDE_Activation = mean(PDE_Activation(:,:,:,:,times_to_plot(its)),4);
    av_APP_Activation = mean(APP_Activation(:,:,:,:,times_to_plot(its)),4);
    
    av_ABM_FullActivation = mean(ABM_FullActivation(:,:,:,:,times_to_plot(its)),4);
    av_PDE_FullActivation = mean(PDE_FullActivation(:,:,:,:,times_to_plot(its)),4);
    av_APP_FullActivation = mean(APP_FullActivation(:,:,:,:,times_to_plot(its)),4);

    mean_ABM_Activation = mean(ABM_Heterogeneity_mean(:,:,:,:,times_to_plot(its)),4);
    var_ABM_Activation = mean(ABM_Heterogeneity_var(:,:,:,:,times_to_plot(its)),4);
    

    mean_TCellAroundDC = mean(ABM_TCellAroundDC_mean(:,:,:,:,times_to_plot(its)),4);
    var_TCellAroundDC = mean(ABM_TCellAroundDC_var(:,:,:,:,times_to_plot(its)),4);
    
    
    C_levels = [0:0.1:1];
        
    
    Cells_Pre_Cluster = ModelParams.NumDCs./NumberOfClusters;

    iter = 0;
    
    
    for pD = 1:length(P_D)
    for pA = 1:length(P_A)
        iter = iter + 1;

        marker_i = marker_style(1);
        
        color_i = its;
        subplot(length(P_D),length(P_A),iter)
        
        hold on;
        plot(log2(Cells_Pre_Cluster) , squeeze(av_ABM_Activation(pA,pD,:)), "-", "Color", cmap(color_i,:), 'LineWidth',2)
        plot(log2(Cells_Pre_Cluster) , squeeze(av_ABM_Activation(pA,pD,:)), marker_i, "Color", cmap(color_i,:), 'LineWidth',1)
        % plot(log2(NumberOfClusters) , squeeze(av_PDE_Activation(pA,pD,:)), 'o-', "Color", yellow_colour, 'LineWidth',2)
        plot(log2(Cells_Pre_Cluster) , squeeze(av_APP_Activation(pA,pD,:)), '--', "Color", font_colour, 'LineWidth',1)
        ax = gca;
        % ax.YScale = 'log';
        hold off;
        axis([min(log2(Cells_Pre_Cluster))-0.5 max(log2(Cells_Pre_Cluster))+0.5 0 1])
        if (iter < 4)
            title(['$\mu_+$ = ', num2str(P_A(pA))], 'Interpreter','latex', 'FontSize',14)
        end
        if(iter > 0)
            xlabel('$\log_2$(Cluster Size)', 'Interpreter','latex', 'FontSize',14)
        end
        if(ismember(iter,[1 4 7]))
           ylabel(['Prop. activated T cells'], 'Interpreter','latex', 'FontSize',14)

           % text(gca, -0.4, 0.5, ['$\mu_-$ = ', num2str(P_D(pD))], ...
           %  'Units', 'normalized', ...
           %  'HorizontalAlignment', 'center', ...
           %  'VerticalAlignment', 'middle', ...
           %  'FontSize',14, 'Interpreter','latex');
        end
        if(ismember(iter,[3 6 9]))
             text(gca, 1.1, 0.5, ['$\mu_-$ = ', num2str(P_D(pD))], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Rotation',90,...
             'FontSize',14, 'Interpreter','latex','Color',font_colour);
         end
    end
    end
    % subplot(length(P_D),length(P_A),1)
    % hold on;
    % ylabel('Proprtion of activated T cells', 'Interpreter','latex', 'FontSize',14)
    % 

end

sgtitle(['T cell activation profile'], 'Interpreter','latex', 'FontSize',16,'Color',font_colour);


formatSpec = '%1.0f';
    
subplot(length(P_D),length(P_A),1)
hold on;
p = [];
str = [];
for its = 1:length(times_to_plot)
    time = times_to_plot(its)*ModelParams.T_final/NumberOfTimeSamples * (1/(60));
    p = [p; plot(-1 , -1, 'o', "Color", cmap(its,:), 'LineWidth',2)];
    str = [str; num2str(time,formatSpec) + " hrs"];
end
p = [p; plot(-1 , -1, '--', "Color", back_colour, 'LineWidth',1)];
str = [str; "Approx."];
legend(p, str,'Location','northwest', 'Interpreter','latex', 'FontSize',12, 'NumColumns',2)


f = gcf;
f.Color = [1 1 1];
f.Position = [251 216 981 650];
% export_fig MainResultsFig_TCellPrifile.png -m2.5

%%
for ii=1:9
    ax = subplot(3,3,ii)
    exportgraphics(ax,['Figs/MainResultsFig_TCellPrifile_subplot_',num2str(ii),'.png'],'Resolution',500)
end
%%

figure;

% times_to_plot = [2:4:NumberOfTimeSamples];
marker_style = ["o-", "square:", "^:"];
cmap = myColour3Gradient(length(times_to_plot), yellow_colour, red_colour, purple_colour);


for its = 1:length(times_to_plot)
    
    
   
    av_ABM_Activation = mean(ABM_Activation(:,:,:,:,times_to_plot(its)),4);
    av_PDE_Activation = mean(PDE_Activation(:,:,:,:,times_to_plot(its)),4);
    av_APP_Activation = mean(APP_Activation(:,:,:,:,times_to_plot(its)),4);
    
    av_ABM_FullActivation = mean(ABM_FullActivation(:,:,:,:,times_to_plot(its)),4);
    av_PDE_FullActivation = mean(PDE_FullActivation(:,:,:,:,times_to_plot(its)),4);
    av_APP_FullActivation = mean(APP_FullActivation(:,:,:,:,times_to_plot(its)),4);

    mean_ABM_Activation = mean(ABM_Heterogeneity_mean(:,:,:,:,times_to_plot(its)),4);
    var_ABM_Activation = mean(ABM_Heterogeneity_var(:,:,:,:,times_to_plot(its)),4);
    

    mean_TCellAroundDC = mean(ABM_TCellAroundDC_mean(:,:,:,:,times_to_plot(its)),4);
    var_TCellAroundDC = mean(ABM_TCellAroundDC_var(:,:,:,:,times_to_plot(its)),4);
    
    
    C_levels = [0:0.1:1];
        
    
    Cells_Pre_Cluster = ModelParams.NumDCs./NumberOfClusters;

    iter = 0;
    
    
    for pD = 1:length(P_D)
    for pA = 1:length(P_A)
        iter = iter + 1;

        marker_i = marker_style(1);
        
        color_i = its;
        subplot(length(P_D),length(P_A),iter)
        
        hold on;
        % plot(log2(Cells_Pre_Cluster) , squeeze(av_ABM_Activation(pA,pD,:)), ":", "Color", cmap(color_i,:), 'LineWidth',2)
        % plot(log2(Cells_Pre_Cluster) , squeeze(av_ABM_Activation(pA,pD,:)), marker_i, "Color", cmap(color_i,:), 'LineWidth',1)
        
        plot(log2(Cells_Pre_Cluster) , squeeze(mean_ABM_Activation(pA,pD,:)), ':', "Color", cmap(color_i,:), 'LineWidth',2)
        plot(log2(Cells_Pre_Cluster) , squeeze(mean_ABM_Activation(pA,pD,:)), marker_i, "Color", cmap(color_i,:), 'LineWidth',1)
        % plot(log2(Cells_Pre_Cluster) , squeeze(mean_ABM_Activation(pA,pD,:)) - squeeze(var_ABM_Activation(pA,pD,:)), '--', "Color", cmap(color_i,:), 'LineWidth',1)
        % plot(log2(Cells_Pre_Cluster) , squeeze(mean_ABM_Activation(pA,pD,:)) + squeeze(var_ABM_Activation(pA,pD,:)), '--', "Color", cmap(color_i,:), 'LineWidth',1)

        hold off;
        % ax.YScale = 'log';

        axis([min(log2(Cells_Pre_Cluster))-0.5 max(log2(Cells_Pre_Cluster))+0.5 0 ModelParams.NumDCs])
        if (iter < 4)
            title(['$\mu_+$ = ', num2str(P_A(pA))], 'Interpreter','latex', 'FontSize',14)
        end
        if(iter > 6)
            xlabel('$\log_2$(Cluster Size)', 'Interpreter','latex', 'FontSize',14)
        end
         if(ismember(iter,[1 4 7]))
            ylabel(['Prop. activated T cells'], 'Interpreter','latex', 'FontSize',14)
         end
         if(ismember(iter,[3 6 9]))
             text(gca, 1.1, 0.5, ['$\mu_-$ = ', num2str(P_D(pD))], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Rotation',90,...
             'FontSize',14, 'Interpreter','latex');
         end
    end
    end
    % subplot(length(P_D),length(P_A),1)
    % hold on;
    % ylabel('Proprtion of activated T cells', 'Interpreter','latex', 'FontSize',14)
    % 

end

sgtitle(['Number of unique DCs a T cell interacts with'], 'Interpreter','latex', 'FontSize',16)


formatSpec = '%1.0f';
    
subplot(length(P_D),length(P_A),1)
hold on;
p = [];
str = [];
for its = 1:length(times_to_plot)
    time = times_to_plot(its)*ModelParams.T_final/NumberOfTimeSamples * (1/(60));
    p = [p; plot(-1 , -1, 'o', "Color", cmap(its,:), 'LineWidth',2)];
    str = [str; num2str(time,formatSpec) + " hrs"];
end
legend(p, str,'Location','northwest', 'Interpreter','latex', 'FontSize',12, 'NumColumns',2)



%%



figure;

cmap2 = myColour3Gradient(length(Cells_Pre_Cluster), yellow_colour , red_colour,  purple_colour );

for its = length(times_to_plot)
  
    Cells_Pre_Cluster = ModelParams.NumDCs./NumberOfClusters;
    
    iter = 0;
    for pD = 1:length(P_D)
    for pA = 1:length(P_A)
        
        iter = iter + 1;
        av_DC_dist = squeeze(mean(ABM_NumbDCAct(pA,pD,:,:,:,times_to_plot(its)),4));


        subplot(length(P_A),length(P_D),iter)
        hold on;
        for jj=[0:32:128]
            plot3(linspace(0,7, 10), jj*ones(1,10),0*ones(1,10), ':' ,'Color',back_colour ,'LineWidth',1.5)
        end

        max_scale = 0; min_scale = 1;
        for jj=1:length(Cells_Pre_Cluster)
            % min_DC = min(find(av_DC_dist(jj,:) ~= 0));
            % if min_DC > 1
            %     min_DC = max_DC - min_DC - 1;
            % end
            % max_DC = max(find(av_DC_dist(jj,:) ~= 0));
            % if max_DC < ModelParams.NumDCs
            %     max_DC = max_DC + 1;
            % end
            % 
            % x_fil = log2(Cells_Pre_Cluster(jj))*ones(1,(max_DC - min_DC)+1);
            % y_fil = min_DC:max_DC;
            % z_fil = [av_DC_dist(jj,min_DC:max_DC)];

            x_fil = log2(Cells_Pre_Cluster(jj))*ones(1,ModelParams.NumDCs+1);
            y_fil = 0:ModelParams.NumDCs;
            z_fil = [av_DC_dist(jj,:)];
            fill3([x_fil fliplr(x_fil)], [y_fil zeros(size(y_fil))], [z_fil z_fil], ...
            cmap2(jj,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none');

            plot3(x_fil, y_fil, z_fil, '-', 'LineWidth',2, 'Color', cmap2(jj,:))

            max_scale = max([max_scale z_fil]);
            min_scale = min([min_scale z_fil(z_fil > 0)]);
        end


        axis([min(log2(Cells_Pre_Cluster)) max(log2(Cells_Pre_Cluster)) 0 ModelParams.NumDCs 0 1.01*max_scale])
        % view(95,30)
        view(80,25)
        zticks([0:0.05:0.2])
        xticks([])
        yticks([0:32:128])

        hold off;
        if (iter < 4)
            title(['$\mu_+$ = ', num2str(P_A(pA))], 'Interpreter','latex', 'FontSize',14)
        end
        if(iter > 6)
            ylabel(['\# DCs to activate'],  'Interpreter','latex', 'FontSize',14)
        end
         if(ismember(iter,[1 4 7]))
            zlabel(['Prop. activated T cells'], 'Interpreter','latex', 'FontSize',14)
         end
         if(ismember(iter,[3 6 9]))
             text(gca, 1.1, 0.5, ['$\mu_-$ = ', num2str(P_D(pD))], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Rotation',90,...
             'FontSize',14, 'Interpreter','latex','Color',back_colour);
         end
    
    end
    end


    % subplot(2,length(P_A),length(P_A)+1)
    % hold on;
    % ylabel('\# DCs to activate', 'Interpreter','latex', 'FontSize',14)
    
    
    formatSpec = '%i';
    
    subplot(length(P_A),length(P_D),1)
    hold on;
    p = [];
    str = [];
    for jj=length(Cells_Pre_Cluster):-1:1
        if jj==0;
            % p = [p; scatter(-1,-1,'Marker','none')];
            % str = [str; "$Log_{2}$(Size)" ];
        else
            p = [p; plot(-1 , -1, '-', "Color", cmap2(jj,:), 'LineWidth',2)];
            str = [str; ""+num2str((Cells_Pre_Cluster(jj)),formatSpec)];
        end
        
    end
    lgd = legend(p, str,'Location','northwest', 'Interpreter','latex', 'FontSize',12, NumColumns=2);
    lgd.ItemTokenSize = [10,5];
    lgd.Title.String = 'Cluster Size';
    lgd.Position = [0.2561 0.7224 0.0721 0.1458];

    % sgtitle('Heterogeneity of T cell activation', 'Interpreter','latex', 'FontSize',16,'Color',back_colour)
end

f = gcf;
f.Color = [1 1 1];
f.Position = [251 216 981 650];
% export_fig MainResultsFig_TCellHeterogeneity.png -m2.5

%%
for ii=1:9
    ax = subplot(3,3,ii)
    exportgraphics(ax,['Figs/MainResultsFig_TCellHeterogeneity_subplot_',num2str(ii),'.png'],'Resolution',500)
end
