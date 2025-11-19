% clc;
% close all;
clear all;

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];

%%
% close all;
% profile on;

ModelParams = struct();

domainBoundary = struct();
domainBoundary.x_max =  9;
domainBoundary.y_max =  9;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 500;
ModelParams.p_move = 50 * 1/16; %(unit conversion)
ModelParams.p_move = 0.5; %(unit conversion)

% ChemoTaxis params
ModelParams.C_sens = 2.0*ModelParams.p_move;             % Taxis sensitivity coefficient
ModelParams.C_sens = 1.0;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 1;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 1;

% T-cells
ModelParams.P_A = 0.5; % Antigen gain rate
ModelParams.P_D = 0.5;  % Antigen loss rate
ModelParams.activatedAge  = 50;

dx = 0.25;
dy = dx;
dt = 0.01;
da = ModelParams.activatedAge/50;

ModelParams.dx_ABM = dx;
ModelParams.dy_ABM = dy;
ModelParams.dt_ABM = dt;            % Time step
ModelParams.da_ABM = da;            % Time step
ModelParams.NumbAgents = 1e4;

ModelParams.dx_PDE = dx;
ModelParams.dy_PDE = dy;
ModelParams.dt_PDE = dt;
ModelParams.da_PDE = da;

ModelParams.plot_traj = true;
recordRate = 31.25/dt;
recordRate = 25/dt;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;


% rndSeed = 6; % for 1 DC and domain -> 9
% rndSeed = 4; % for 1 DC and domain -> 20
% rndSeed = 5; % for 2 DC and domain -> 20
% rndSeed = 4; % for 4 DC and domain -> 20
% rndSeed = 4; % for 7 DC and domain -> 20
rndSeed = 3

[(ModelParams.p_move)*(4*dt)/(dx.^2) (ModelParams.C_sens)*(4*dt)/(dx.^2) (ModelParams.P_A/da)*dt (ModelParams.P_A/ModelParams.activatedAge)*dt (ModelParams.P_D/da)*dt ]


[u, C_pde, dCdx, dCdy, Ic, A_Ic, Inds_Ic_st, params, A] = Phenotype_PDE_SetUp(rndSeed,domainBoundary,ModelParams);
activation_proportion_pheno = zeros(params.nt/recordRate,1);
Dx = params.D;
chix = params.C_chi;
rho_plus = params.P_A/params.activatedAge;
rho_minus = params.P_D/params.activatedAge;

inv_dx = 1.0/params.dx;
inv_dy = 1.0/params.dy;
inv_da = 1.0/params.da;
dt = params.dt;

totalsize = params.Nx*params.Ny*params.Na;

i_w = true(totalsize,1); i_w(Inds_Ic_st.inds_west)  = false;
i_e = true(totalsize,1); i_e(Inds_Ic_st.inds_east)  = false;
i_s = true(totalsize,1); i_s(Inds_Ic_st.inds_south) = false;
i_n = true(totalsize,1); i_n(Inds_Ic_st.inds_north) = false;
i_c = true(totalsize,1); i_c(Inds_Ic_st.inds)       = false;

dims = int32([params.Nx, params.Ny, params.Na]);

x = linspace(0, params.Lx, params.Nx);
y = linspace(0, params.Ly, params.Ny);
a = linspace(0, params   .La, params.Na);


[walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
activatedCells_ABM = zeros(params_ABM.nt/recordRate,1);
seed = rndSeed;

figure;
% subplot(1,2,1)
imagesc((DC_model.BoundaryDC))
% imagesc((C) - C_pde(:,:,1)')
%%


% compute analytic values for the proportion of activation and mean ac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, C_an, Ic_an, A_Ic_an, params_an, Inds_Ic_an] = Analytic_SetUp(rndSeed,domainBoundary,ModelParams);
a_an = linspace(da, params_an.La-da, params_an.Na-2);
a_an = linspace(0, params_an.La, params_an.Na);

u_SS_sol = exp(params_an.C_chi/params_an.D .* C_an);

u_SS_sol(sub2ind([params_an.Ny, params_an.Nx], Inds_Ic_an(:,1), Inds_Ic_an(:,2))) = 0;
u_SS_sol = u_SS_sol/sum(u_SS_sol(:));

a_SS_sol = A_Ic_an(:,:,end) .* u_SS_sol;

Prob_A = sum(a_SS_sol(:));
Prob_Ac = 1 - Prob_A;


LengthScale_A = sqrt((nnz(A_Ic_an)/(ModelParams.numberOfClusters))*ModelParams.dx_ABM^2/pi);
if ismember(ModelParams.numberOfClusters, [2].^[0 2 4 6 8])
    if ModelParams.numberOfClusters == 1
        y_mx_i = 0.5*domainBoundary.y_max;
        x_mx_i = 0.5*domainBoundary.x_max;
        LengthScale_Ac = (1/2)*( (y_mx_i-LengthScale_A) + sqrt(y_mx_i^2+x_mx_i^2)-LengthScale_A);
    else
        y_mx_i = 0.5*domainBoundary.y_max/sqrt(ModelParams.numberOfClusters);
        x_mx_i = 0.5*domainBoundary.x_max/sqrt(ModelParams.numberOfClusters);
        LengthScale_Ac = (1/2)*( (y_mx_i-LengthScale_A) + sqrt(y_mx_i^2+x_mx_i^2)-LengthScale_A);
    end
else
    y_mx_i = 0.5*domainBoundary.y_max/sqrt(ModelParams.numberOfClusters/2);
    x_mx_i = 0.5*domainBoundary.x_max/sqrt(ModelParams.numberOfClusters*2);
    lngth = max(y_mx_i,x_mx_i);
    LengthScale_Ac = (1/2)*( (lngth-LengthScale_A) + sqrt(y_mx_i^2+x_mx_i^2)-LengthScale_A);
end
sigma_val = 6*(ModelParams.p_move) /( (LengthScale_Ac*LengthScale_A));

chi_sense = ModelParams.C_sens;
D_move    = ModelParams.p_move;
gain_rate = ModelParams.P_A ;
loss_rate = ModelParams.P_D ;
max_A     = ModelParams.activatedAge;

alpha_val = sigma_val *  Prob_A  * max_A / loss_rate; 
beta_val =  sigma_val *  Prob_Ac * max_A / gain_rate;
            
U_scale = (1./beta(alpha_val,beta_val)) * (1/(max_A).^(alpha_val+beta_val)) .* (alpha_val+beta_val)./(alpha_val+ (gain_rate/loss_rate) * beta_val);
U_A =                                     a_an .*                           a_an.^(alpha_val-1) .*(max_A - a_an).^(beta_val-1);
U_Ac =  (gain_rate/loss_rate) .* (max_A - a_an) .* a_an.^(alpha_val-1) .*(max_A - a_an).^(beta_val-1);
kappa_1 = alpha_val;
kappa_2 = beta_val;

mean_act = (max_A * kappa_1)/(kappa_1+kappa_2+1) * ((1/(kappa_1 + (gain_rate/loss_rate)*kappa_2)) + 1);
mean_act_level = mean_act/max_A;
u_of_a = U_scale.*(U_A+U_Ac)*da;

activation_proportion = mean_act_level;
expected_activation_level = mean_act;

disp(["kappa_1 = ", num2str(kappa_1), " kappa_2 = ", num2str(kappa_2)])

subplot(1,4,3)
hold on;
plot(a_an,u_of_a, LineWidth=2.5)



%%
pause(1)
cmap = myColour3Gradient(params.nt/recordRate + 1,yellow_colour, red_colour, purple_colour);
c_it = 0;
tic;

figure;

for n=0:(params_ABM.nt/recordRate)
    seed = seed + 1;

    if true
        % clf;
        c_it = c_it + 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        u_slice = squeeze(sum(u,3));
        % surf(x, y, u_slice', 'EdgeColor','none');
        % colorbar;
        % view(2)
        % ax1.XTick = [];
        % ax1.YTick = [];
        % xlabel('x', 'FontSize',14, 'Interpreter','latex');
        % ylabel('y', 'FontSize',14, 'Interpreter','latex');
        % title('PS-PDE T cell density', 'FontSize',14, 'Interpreter','latex');
        % axis([0 params.Lx 0 params.Ly])

        av_a = zeros(params.Nx, params.Ny);
        for ix=1:params.Nx
            for iy=1:params.Ny
                av_a(ix,iy) = a*squeeze(u(ix,iy,:));
            end
        end
        av_a = av_a/(params.activatedAge*sum(u(:)));
        % av_a = av_a/(sum(u(:)));

        % 
        % ax2 = subplot(2,3,2);
        % surf(x, y, av_a', 'EdgeColor','none');
        % colorbar;
        % view(2)
        % ax2.XTick = [];
        % ax2.YTick = [];
        % axis([0 params.Lx 0 params.Ly])
        % xlabel('x', 'FontSize',14, 'Interpreter','latex');
        % ylabel('y', 'FontSize',14, 'Interpreter','latex');
        % title('PS-PDE antigen density', 'FontSize',14, 'Interpreter','latex');        
        

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        

        A_abm = zeros(params_ABM.Ny, params_ABM.Nx);
        U_abm = zeros(params_ABM.Ny, params_ABM.Nx);

        for i = 1:params_ABM.num_walkers
            id_x_i = round(walker_positions(i,1)/params_ABM.dx) + 1;
            id_y_i = round(walker_positions(i,2)/params_ABM.dy) + 1;
        
            U_abm(id_y_i,id_x_i) = U_abm(id_y_i,id_x_i) + 1;
            A_abm(id_y_i,id_x_i) = A_abm(id_y_i,id_x_i) + walker_activation(i);
        end
        A_plot = zeros(size(U_abm));
        for ii=1:size(U_abm,1)
            for jj=1:size(U_abm,2)
                if U_abm(jj,ii) > 0 
                    A_plot(jj,ii) = A_abm(jj,ii)/(params_ABM.activatedAge*U_abm(jj,ii));
                end
            end
        end
            
        Ic = DC_model.Ic;
        U_plot = (U_abm/params_ABM.num_walkers);
        U_plot_1 = U_plot;
        tmp_act = (A_plot).*U_plot;
        
        % for x_max = 10
        % x_cont_1 = 0.002;
        % x_cont_2 = 0.004;
        
        % for x_max = 20
        % x_cont_1 = 0.95e-3;   % dashed
        % x_cont_2 = 6e-4; % dotted
        x_cont_1 = 1e-3;   % dashed
        x_cont_2 = 5e-4; % dotted


        % x_cont_1 = 0.0022;
        % x_cont_2 = 0.003;
        
        
        % subplot(3,4,1)
        % surf(x,y,u_slice','edgecolor','none')
        % view(2)
        % subplot(3,4,5)
        % surf(x,y,U_plot_1,'edgecolor','none')
        % view(2)
        % 
        % subplot(3,4,9)
        % mean_abm_activation = mean(walker_activation);
        % NumberDCsToActivate_tmp = zeros(ModelParams.NumDCs,1);
        % counts = 0;
        % dcs_act = [];
        % for ii=1:size(DCLingerTime,1)
        %     if walker_activation(ii) >= 0
        %         ii_NumberDCsToActivate = nnz(DCLingerTime(ii,:));
        %         if ii_NumberDCsToActivate > 0
        %             dcs_act = [dcs_act ii_NumberDCsToActivate];
        %             counts = counts + 1;
        %             NumberDCsToActivate_tmp(ii_NumberDCsToActivate) = NumberDCsToActivate_tmp(ii_NumberDCsToActivate) + 1;
        %         end
        %     end
        % end
        % NumberDCsToActivate = NumberDCsToActivate_tmp(:)/counts;
        % plot(NumberDCsToActivate,'LineWidth',2)

        ax1 = subplot(1,4,1);
        hold on;
        contour(ax1, x,y,U_plot_1, [x_cont_1, x_cont_1], 'linewidth',2, 'LineStyle', '-', 'EdgeAlpha', 0.8 , 'EdgeColor', cmap(c_it,:));
        contour(ax1, x,y,u_slice', [x_cont_1, x_cont_1], 'linewidth',2, 'LineStyle', '--', 'EdgeAlpha', 0.8, 'EdgeColor', back_colour);

        contour(ax1, x,y,U_plot_1, [x_cont_2, x_cont_2], 'linewidth',2, 'LineStyle', ':', 'EdgeAlpha', 0.8 , 'EdgeColor', cmap(c_it,:));
        contour(ax1, x,y,u_slice', [x_cont_2, x_cont_2], 'linewidth',2, 'LineStyle', ':', 'EdgeAlpha', 0.8, 'EdgeColor', back_colour);

        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        view(0,90)
        title('T cell density', 'FontSize',16 ,'Interpreter','latex')
        xlabel('x', 'FontSize',14 ,'Interpreter','latex')
        ylabel('y', 'FontSize',14 ,'Interpreter','latex')

        
        % a_cont = 1e-4;
        % ax2 = subplot(1,4,2);
        % hold on;
        % contour(ax2, x,y,tmp_act, [a_cont, a_cont], 'linewidth',2, 'LineStyle', '-', 'EdgeAlpha', 0.8, 'EdgeColor', cmap(c_it,:));
        % contour(ax2, x,y,av_a',   [a_cont, a_cont], 'linewidth',2, 'LineStyle', '--', 'EdgeAlpha', 0.8, 'EdgeColor', back_colour);
        % axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        % view(0,90)
        % ax2.XTick = [];
        % ax2.YTick = [];
        


        % ax3 = subplot(1,4,4);
        % hold on
        % plot(ax3, recordRate*dt*(0:n), [0; activatedCells_ABM(1:n)],'-','LineWidth',2.5, 'Color', purple_colour);
        % plot(ax3, recordRate*dt*(0:n), [0; activation_proportion_pheno(1:n)], '--', 'LineWidth',2.5, 'Color',yellow_colour)
        % 
        % plot(ax3, 0:recordRate, activation_proportion*ones(recordRate+1,1), ':' ,'LineWidth', 3, 'Color', blue_colour )
        % 
        % xlabel('Time', 'FontSize',14, 'Interpreter','latex'); 
        % axis([0 params.nt *dt 0 1])
        % % title(['Proportion of activation ABM: ', num2str(activatedCells_ABM((n+1))), ', PDE: ' ,num2str(activation_proportion_pheno((n+1)))], 'FontSize',14, 'Interpreter','latex')
        % title('Proportion of T cell activation', 'FontSize',16 ,'Interpreter','latex')
        % 
        % 
        [N_A, A_edges] = histcounts(walker_activation,-(da*0.5):da:(params_ABM.activatedAge+0.5*da));
        PDE_dist =  squeeze(sum(u,[1,2]));
        
        % ax4 = subplot(1,4,4);
        % hold on;
        % plot(ax4, a, PDE_dist, 'color', back_colour , 'LineWidth',1.5)
        % plot(ax4, A_edges(2:end)-(da*0.5), N_A/(ModelParams.NumbAgents),'--','LineWidth',2, 'Color', cmap(c_it,:))
        % % plot(ax4, mean(walker_activation)*[1 1], [0 1.25*max(N_A)], '--','color', blue_colour ,'LineWidth', 2 )
        % % plot(ax4, [expected_activation_level expected_activation_level], [0 1.25*max(N_A)], ':','color', yellow_colour ,'LineWidth', 2 )
        % % plot(ax4, [a*PDE_dist a*PDE_dist], [0 1.25*max(N_A)], '--','color', red_colour ,'LineWidth', 2 )
        % 
        % plot(ax4, a_an, u_of_a, ':' ,'color', yellow_colour , 'LineWidth',3)
        % 
        % axis([0 params_ABM.activatedAge 0 1.1*max([N_A/(ModelParams.NumbAgents)]) ])
        % xlabel('Time', 'FontSize',14, 'Interpreter','latex'); 
        % % axis([0 params_ABM.activatedAge])
        % title(['Antigen distribution ABM: ', num2str((N_A(end)/ModelParams.NumbAgents)), ', PDE: ' , num2str((PDE_dist(end-1)+PDE_dist(end)))], 'FontSize',14, 'Interpreter','latex')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        subplot(1,4,2)
        time_ii = recordRate*dt*n;
        hold on;
        x_fil = time_ii*ones(length(N_A(1:end))+1,1)';
        y_fil = A_edges(1:end)-0.5*da;
        z_fil = [0 N_A/(ModelParams.NumbAgents)];
        % fill3([x_fil fliplr(x_fil)], [y_fil zeros(size(y_fil))], [z_fil z_fil], ...
        % cmap(c_it,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none');

        plot3(time_ii*ones(length(N_A(1:end)),1), A_edges(1:end-1)+0.5*da, N_A/(ModelParams.NumbAgents), '-', 'LineWidth',2, 'Color', cmap(c_it,:))
        plot3((time_ii-0.1)*ones(length(a),1), a ,PDE_dist,                                   '--' ,'LineWidth',1.5, 'Color',back_colour)


        axis([0 params_ABM.nt*dt 0 params_ABM.activatedAge 0 0.08])
        title('Stimulation distribution', 'FontSize',16 ,'Interpreter','latex')
        xlabel('t', 'FontSize',14 ,'Interpreter','latex')
        ylabel('a', 'FontSize',14 ,'Interpreter','latex')
        view(-50,25)

       


        % sgtitle(['time = ', num2str(params.dt*n), ' mass = ', num2str(sum(u(:)))])
        drawnow;
        
    end

    % % % % % u = computePhenotypeModel(u, dCdx, dCdy, A, Inds_Ic_st, A_Ic, params);
     u = computePhenotypeModelTime_mex( ...
        u,    C_pde,    A,    A_Ic,  ...
        params.C_chi, params.D,   1.0/params.dx,   1.0/params.dy,    ...
        params.P_A/params.activatedAge, params.P_D/params.activatedAge, ...
        params.activatedAge, 1.0/params.da,  params.dt, ...
        i_e,  i_w,  i_n,  i_s,   i_c, ...
        dims, recordRate);
    % % % % activation_proportion_pheno(n) = squeeze(sum(u,[1,2]))'*a';
    activation_proportion_pheno(n+1) = squeeze(sum(u,[1,2]))'*a'/params.activatedAge;
    
    % Call the C MEX function
    [walker_positions, walker_activation, C, DCLingerTime] = ...
        computeABMModelTime_mex( ...
            walker_positions, ...       % 1
            walker_activation, ...      % 2
            C, ...                      % 3
            DCLingerTime, ...           % 4
            DC_model.Ic, ...            % 5
            DC_model.BoundaryDC, ...    % 6
            params_ABM.num_walkers, ... % 7
            params_ABM.dx, ...          % 8
            params_ABM.dy, ...          % 9
            params_ABM.Nx, ...          % 10
            params_ABM.Ny, ...          % 11
            params_ABM.p_move, ...      % 12
            params_ABM.C_chi, ...       % 13
            params_ABM.activatedAge, ...% 14
            params_ABM.P_A, ...         % 15
            params_ABM.P_D, ...         % 16
            params_ABM.dt, ...          % 17
            params_ABM.da, ...          % 18
            ModelParams.NumDCs, ...
            recordRate ...
        );

    activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*params_ABM.num_walkers);


end

ax4 = subplot(1,4,3);
hold on;
plot(ax4, A_edges(2:end)-(da*0.5), N_A/(ModelParams.NumbAgents),'-','LineWidth',2.5, 'Color', purple_colour)
plot(ax4, a, PDE_dist, '--', 'color', yellow_colour , 'LineWidth',2.5)
plot(ax4, [expected_activation_level expected_activation_level], [0 1.25*max(N_A)], ':','color', blue_colour ,'LineWidth', 3 )


plot(ax4, mean(walker_activation)*[1 1], [0 1.25*max(N_A)], '-','color', purple_colour ,'LineWidth', 2.5 )
plot(ax4, [a*PDE_dist a*PDE_dist], [0 1.25*max(N_A)], '--','color', yellow_colour ,'LineWidth', 2.5 )
plot(ax4, a_an, u_of_a, ':' ,'color', blue_colour , 'LineWidth',3)

axis([0 params_ABM.activatedAge 0 0.08])
xlabel('a', 'FontSize',16, 'Interpreter','latex'); 
% axis([0 params_ABM.activatedAge])
title('Stimulation distribution (steady state)', 'FontSize',14 ,'Interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






ax3 = subplot(1,4,4);
hold on


plot(ax3, recordRate*dt*(0:n), [0; activatedCells_ABM(1:n)],'-','LineWidth',2.5, 'Color', purple_colour);
plot(ax3, recordRate*dt*(0:n), [0; activation_proportion_pheno(1:n)], '--', 'LineWidth',2.5, 'Color',yellow_colour)

plot(ax3, 0:params_ABM.nt*dt, activation_proportion*ones(length(0:params_ABM.nt*dt),1), ':' ,'LineWidth', 3, 'Color', blue_colour )

scatter(-1,-1,'Marker','none')
scatter(-1,-1,'Marker','none')
% scatter(-1,-1,'Marker','none')

xlabel('Time', 'FontSize',14, 'Interpreter','latex'); 
axis([0 params.nt*dt 0 1])
% title(['Proportion of activation ABM: ', num2str(activatedCells_ABM((n+1))), ', PDE: ' ,num2str(activation_proportion_pheno((n+1)))], 'FontSize',14, 'Interpreter','latex')
title('Proportion of T cell activation', 'FontSize',16 ,'Interpreter','latex')


lg_str = {'ABM', 'PS-PDE','Approx.',...
           %  ['$\kappa_{1} \approx' , num2str(round(kappa_1,2)) ,'$'], ...
           % ['$ \kappa_{2} \approx' ,num2str(round(kappa_2,2)),'$'], ...
           };

% lg_str = {'ABM', 'PS-PDE','Approx.'};

% lg_str = [{['$A_{max}=' , num2str(max_A),'$'] , ['$\kappa_{1} \approx' , num2str(round(kappa_1,2)) ,'$'], ...
%            ['$ \kappa_{2} \approx' ,num2str(round(kappa_2,2)),'$'], ...
%            'ABM', 'PS-PDE','Approx.'}];

% lg_str = [{['  \mu_+  = ' , num2str(gain_rate) ,'  '],['  \mu_-  = ' , num2str(loss_rate),'  '], ...
%            ['A_{max}= ' , num2str(max_A),''] , ['    \chi   = ' , num2str(chi_sense),''] ,...
%            'Analytic', 'ABM', 'PS-PDE'}]
% lg = legend(lg_str, 'FontSize',16 , 'location','southoutside',NumColumns=2);
% lg = legend('Analytic', 'ABM', 'PS-PDE', 'FontSize',16 ,'Interpreter','latex', 'location','southoutside');

lg = legend(lg_str, 'FontSize',16 ,'Interpreter','latex',NumColumns=1);


% ax3 = subplot(1,4,4);
% hold on;
% lg_str = [{'ABM', 'PS-PDE','Approx.'}];
% lg = legend(lg_str, 'FontSize',16 ,'Interpreter','latex',NumColumns=1);
toc;

% profile viewer
% profile off;

%%

% figure;
% 
% u_SS_sol = exp(params.C_chi/params.D .* C)';
% 
% u_SS_sol(i_c(1:(params.Nx*params.Ny))==0) = 0;
% u_SS_sol = u_SS_sol/sum(u_SS_sol(:));
% 
% subplot(2,2,1)
% surf(x, y, u_SS_sol', 'EdgeColor','none');
% colorbar;
% view(2)
% XTick = [];
% YTick = [];
% 
% 
% a_SS_sol = A_Ic(:,:,end) .* u_SS_sol;
% 
% subplot(2,2,2)
% surf(x, y, a_SS_sol', 'EdgeColor','none');
% colorbar;
% view(2)
% XTick = [];
% YTick = [];
% 
% Prob_A = sum(a_SS_sol(:));
% Prob_Ac = 1 - Prob_A;
% 
% activation_proportion = (params.P_A*Prob_A)/( params.P_A*Prob_A + params.P_D*Prob_Ac )
% expected_activation_level = params_ABM.activatedAge * activation_proportion
% 
% 
% sigma_val = 1;
% alpha_val = sigma_val*params_ABM.activatedAge * Prob_A / params.P_A; 
% beta_val = sigma_val*params_ABM.activatedAge * Prob_Ac / params.P_D;
% 
% U_A = a.^beta_val .*(params_ABM.activatedAge - a).^(alpha_val-1);
% subplot(2,2,3)
% plot(a,U_A, LineWidth=2.5)

%%

subplot(1,4,1)
yticks([0:2.5:10])
xticks([0:2.5:10])

subplot(1,4,2)
axis([0 params_ABM.nt*dt 0 params_ABM.activatedAge 0 0.12])
xticks([0:100:500])
yticks([0:10:50])
zticks([0:0.025:0.11])

subplot(1,4,3)
axis([0 params_ABM.activatedAge 0 0.12])
yticks([0:0.025:0.12])
xlabel('a', 'FontSize',16, 'Interpreter','latex'); 



f = gcf;
f.Position = [9 180 1503 313];
f.Color = [1 1 1];
% export_fig 1Cells.png -m2.5

%%
figure;

mean_abm_activation = mean(walker_activation);

NumberDCsToActivate_tmp = zeros(ModelParams.NumDCs,1);
counts = 0;
dcs_act = [];
for ii=1:size(DCLingerTime,1)
    if walker_activation(ii) >= mean_abm_activation
        ii_NumberDCsToActivate = nnz(DCLingerTime(ii,:));
        if ii_NumberDCsToActivate > 0
            dcs_act = [dcs_act ii_NumberDCsToActivate];
            counts = counts + 1;
            NumberDCsToActivate_tmp(ii_NumberDCsToActivate) = NumberDCsToActivate_tmp(ii_NumberDCsToActivate) + 1;
        end
    end
end
NumberDCsToActivate = NumberDCsToActivate_tmp(:)/counts;

plot(NumberDCsToActivate,'LineWidth',2)
title(['ClsuterSize = ', num2str(ModelParams.NumDCs/ModelParams.numberOfClusters)])
                                                