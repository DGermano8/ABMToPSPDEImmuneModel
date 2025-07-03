clear all;
% close all;

addpath(genpath('UserFunctions'))

ModelParams = struct();
rndSeed = 107;

domainBoundary = struct();
domainBoundary.x_max =  7;
domainBoundary.y_max =  7;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 100;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 1;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 0.2; % Antigen gain rate
ModelParams.P_D = 0.8;  % Antigen loss rate
ModelParams.activatedAge  = 100;

ModelParams.dx_ABM = 0.5;
ModelParams.dy_ABM = 0.5;
ModelParams.dt_ABM = 0.05;            % Time step
ModelParams.NumbAgents = 50000;

ModelParams.plot_traj = true;
recordRate = 10;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;

ModelParams.p_move*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2)
ModelParams.C_sens*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2)

%%
    purple_colour = 1/255*[131, 96, 150];
    red_colour = 1/255*[237, 123, 123];
    yellow_colour = 1/255*[240, 184, 110];
    blue_colour = 1/255*[106 123 237];
    green_colour = 1/255*[77 150 138];
    back_colour = 1/255*[56 56 56];

    cmap = myColour3Gradient(1+ModelParams.T_final/(ModelParams.t_plot*ModelParams.dt_ABM),  yellow_colour, red_colour, purple_colour);


    [walker_positions, walker_activation, C, DCLingerTime, DC_model, params] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
        
    if ModelParams.plot_traj
        f = figure;
        f.Position = [106 360 1360 360];
    end

    % h = waitbar(0, 'Processing...'); % Create progress bar

    
    activatedCells = zeros(params.nt+1,1);
    tic;
    for n = 0:params.nt

        [walker_positions, walker_activation, U, C, Ap, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params);
        activatedCells(n+1) = sum(walker_activation(:))/(params.activatedAge*sum(U(:)));
        if mod(n, (1/ModelParams.dt_ABM) ) == 0
            % waitbar(n/params.nt, h, sprintf('Progress: %1.2f%%', 100*n/params.nt)); % Update progress
        end

        if mod(n, ModelParams.t_plot) == 0 && ModelParams.plot_traj
            iter = n/ModelParams.t_plot + 1;

            x = linspace(0, params.Lx, params.Nx);
            y = linspace(0, params.Ly, params.Ny);
            [X, Y] = meshgrid(x, y);

            % clf;
            total_density = sum(U(:));
            title_string = ['Time = ', num2str(n * params.dt), ', total activation = ', num2str(100*activatedCells(n+1),'%.1f'), '%', ', Density = ', num2str(total_density)];
            % sgtitle(title_string,'Color',1/255*[217 217 217],'FontSize',18);
            sgtitle(title_string,'FontSize',18);
            
            A_abm = zeros(size(U));
            U_abm = zeros(size(U));
            for i = 1:params.num_walkers
                id_x_i = walker_positions(i,1)/params.dx + 1;
                id_y_i = walker_positions(i,2)/params.dy + 1;
            
                U_abm(id_y_i,id_x_i) = U_abm(id_y_i,id_x_i) + 1;
                A_abm(id_y_i,id_x_i) = A_abm(id_y_i,id_x_i) + walker_activation(i);
            end
            A_plot = zeros(size(U_abm));
            for ii=1:size(U_abm,1)
                for jj=1:size(U_abm,2)
                    if U_abm(jj,ii) > 0 
                        A_plot(jj,ii) = A_abm(jj,ii)/(params.activatedAge*U_abm(jj,ii));
                    end
                end
            end
            U_plot = (U_abm/params.num_walkers);
            clf;
            % U_plot = (U/(sum(U(:))));
            subplot(1,3,1);
            % hold on;
            surf(X, Y, U_plot,'EdgeAlpha',0);
            % contour(X,Y,U_plot, [5.25*10^-3, 5.25*10^-3], 'LineWidth',1.5, 'EdgeColor',cmap(iter,:))
            xlabel('x'); ylabel('y'); title('U');
            c1 = colorbar;
            axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
            % colormap(ax1,cmap2);
            % caxis([0 max(U_plot(:)) ])
            view(0,90)
            ax1.XTick = [0:10:50];
            ax1.YTick = [0:10:50];
            c1.Label.String = 'T-cell density';
            c1.Label.FontSize = 12;
    
            subplot(1,3,2);
            % hold on;
            % % contour(X,Y,(Ap)/sum(U(:)), [2*10^-4, 2*10^-4], 'LineWidth',1.5, 'EdgeColor',cmap(iter,:))
            % contour(X,Y,A_plot, [2*10^-4, 2*10^-4], 'LineWidth',1.5, 'EdgeColor',cmap(iter,:))
            surf(X, Y, A_plot,'EdgeAlpha',0);
            xlabel('x'); ylabel('y'); title('A');
            c2 = colorbar; % caxis([0 1]);
            axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
            % colormap(ax2,cmap2);
            % caxis([0 max(max( (Ap(:)))/sum(U(:)) , 10^-6) ])
            view(0,90)
            ax2.XTick = [0:10:50];
            ax2.YTick = [0:10:50];
            c2.Label.String = 'T-cell activation';
            c2.Label.FontSize = 12;

            subplot(1,3,3);
            plot( 0:n, activatedCells(1:(n+1)), 'Color', blue_colour, ...
                'LineWidth',2.5);
            xlabel('Time'); ylabel('Total activation');
            axis([0 params.nt 0 1])

            % MakeDark();
            drawnow;

        end
    end
    close(h);

    for jj=1:length(walker_activation)
        if (walker_activation(jj) < params.activatedAge)
            DCLingerTime(jj,:) = 0;
        end
    end
   
    
toc;

%%

% 
% purple_colour = 1/255*[131, 96, 150];
% red_colour = 1/255*[237, 123, 123];
% yellow_colour = 1/255*[240, 184, 110];
% blue_colour = 1/255*[106 123 237];
% green_colour = 1/255*[77 150 138];
% back_colour = 1/255*[56 56 56];
% 
% figure;
% plot(0:ModelParams.dt_ABM:ModelParams.T_final, activatedCells,'LineWidth',3.5)
% axis([0 ModelParams.T_final 0 1])
% 
% 
% NumberDCsToActivate = zeros(ModelParams.NumDCs,1);
% for ii=1:size(DCLingerTime,1)
%     if walker_activation(ii) == ModelParams.activatedAge
%         ii_NumberDCsToActivate = nnz(DCLingerTime(ii,:));
%         if ii_NumberDCsToActivate > 0
%             NumberDCsToActivate(ii_NumberDCsToActivate) = NumberDCsToActivate(ii_NumberDCsToActivate) + 1;
%         end
%     end
% end
% 
% 
% 
% 
% 
% f = figure;
% f.Position = [476 360 358 259];
% % Plot the histogram using bar
% if( ModelParams.P_A  <= 0.25)
%     barcolor = red_colour;
% elseif ( ModelParams.P_A  < 0.75 && ModelParams.P_A > 0.25)
%     barcolor = purple_colour;
% else
%     barcolor = blue_colour;
% end
% bar(1:ModelParams.NumDCs, NumberDCsToActivate, 'FaceColor', barcolor, 'EdgeColor', [0.8 0.8 0.8] );
% 
% % Add labels and title
% xlabel('# unique DC interactions');
% % ylabel('Counts');
% mean_DCs = sum(NumberDCsToActivate.*(1:length(NumberDCsToActivate))')/(sum(NumberDCsToActivate));
% median_Count = round(sum(NumberDCsToActivate)*0.5);
% 
% median_val = 0;
% iter = 0;
% for ii=1:length(NumberDCsToActivate)
%     iter = iter + NumberDCsToActivate(ii);
%     if iter > median_Count;
%         median_val = ii
%         break;
%     end
% end
% 
% MakeDark();




%%

