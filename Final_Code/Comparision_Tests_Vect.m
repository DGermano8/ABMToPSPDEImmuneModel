clear all;
% close all;

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];


addpath(genpath('UserFunctions'))

%%

ModelParams = struct();
rndSeed = 107;

domainBoundary = struct();
domainBoundary.x_max =  7;
domainBoundary.y_max =  7;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 2000;
ModelParams.p_move = 1.0;

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
ModelParams.activatedAge  = 50;

% ModelParams.dx_ABM = 0.25;
% ModelParams.dy_ABM = 0.25;
% ModelParams.dt_ABM = 0.01;            % Time step
% ModelParams.NumbAgents = 10000;
% 
% ModelParams.dx_PDE = 0.25;
% ModelParams.dy_PDE = 0.25;
% ModelParams.dt_PDE = 0.01;

ModelParams.dx_ABM = 1/3;
ModelParams.dy_ABM = 1/3;
ModelParams.dt_ABM = 0.025;   % Time step
ModelParams.NumbAgents = 2000;

ModelParams.dx_PDE = 1/3;
ModelParams.dy_PDE = 1/3;
ModelParams.dt_PDE = 0.025;

ModelParams.plot_traj = true;
recordRate = 25;
ModelParams.t_plot = recordRate/ModelParams.dt_ABM;

if (ModelParams.p_move*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    ModelParams.C_sens*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    (4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) == 1)

    disp('Choose smaller dx and dt for ABM')
else
    disp('Nice params!')
end

h = waitbar(0, 'Processing...'); % Create progress bar

%%
x1 = 0:ModelParams.dt_PDE:ModelParams.T_final'; % First x-mesh
x2 = 0:ModelParams.dt_ABM:ModelParams.T_final'; % Second x-mesh


% P_A_vals =  [0.01:0.02:0.09 0.1:0.05:0.35 0.4:0.1:0.9];
% P_D_vals =  [0.01:0.02:0.09 0.1:0.05:0.35 0.4:0.1:0.9];
% C_sens_vals = [0:0.25:1];

P_A_vals =  linspace(0,1,11);
P_D_vals =  linspace(0,1,11);
C_sens_vals = [0.5];

ABM_Activation = zeros(length(P_A_vals), length(P_D_vals), length(C_sens_vals));
PDE_Activation = zeros(length(P_A_vals), length(P_D_vals), length(C_sens_vals));

ABM_FullActivation = zeros(length(P_A_vals), length(P_D_vals), length(C_sens_vals));
PDE_FullActivation = zeros(length(P_A_vals), length(P_D_vals), length(C_sens_vals));

ABM_Activation_Time = zeros(length(P_A_vals), length(P_D_vals), length(C_sens_vals),length(x2));
PDE_Activation_Time = zeros(length(P_A_vals), length(P_D_vals), length(C_sens_vals),length(x1));


error_measure = zeros(length(P_A_vals), length(P_D_vals), length(C_sens_vals));

iter=0;
tic;
for kk=1:length(C_sens_vals)
    ModelParams.C_sens = C_sens_vals(kk);

    for jj=1:length(P_D_vals)
        ModelParams.P_D = P_D_vals(jj);

        for ii=1:length(P_A_vals)
            ModelParams.P_A = P_A_vals(ii);

            iter = iter + 1;
            prog = (iter)/(length(P_D_vals)*length(P_A_vals)*length(C_sens_vals));
            
            % [walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
            % activatedCells_ABM = zeros(params_ABM.nt+1,1);
            % for n = 0:params_ABM.nt
            %     [walker_positions, walker_activation, U, C, Ap, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
            %     activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*sum(U(:)));
            %     % activatedCells_ABM(n+1) = sum(( walker_activation == ModelParams.activatedAge))/(length(walker_activation));
            % end
            % ABM_Activation_Time(ii,jj,kk,:) = activatedCells_ABM;
            % ABM_Activation(ii,jj,kk) = activatedCells_ABM(end);
            % ABM_FullActivation(ii,jj,kk) = nnz(walker_activation == params_ABM.activatedAge)/length(walker_activation);

            [U, A,  U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params_PDE] = PDE_SetUp(rndSeed,domainBoundary,ModelParams);
            activatedCells_PDE = zeros(params_PDE.nt+1,1);
            for n = 0:params_PDE.nt
                [U, A, U_a, A_a] = computePDEModel(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy, params_PDE);
                tmp = (A+A_a);
                % tmp = (A+A_a)/params_PDE.activatedAge;
                % tmp = (A_a)/params_PDE.activatedAge;
                activatedCells_PDE(n+1) = ( sum( tmp(:) ) );
            end
            PDE_Activation_Time(ii,jj,kk,:) = activatedCells_PDE/ModelParams.activatedAge;
            PDE_Activation(ii,jj,kk) = activatedCells_PDE(end)/ModelParams.activatedAge;
            PDE_FullActivation(ii,jj,kk) = sum(A_a(:))/ModelParams.activatedAge;
            
            y1 = activatedCells_PDE'; % First curve
            % y2 = activatedCells_ABM'; % Second curve
            
            % area_between = (computeL1NormError(x1,y1,x2,y2)/ModelParams.T_final)*100;
            % error_measure(ii,jj,kk) = area_between;

            % L2Norm = computeL2NormError(x1,y1,x2,y2);
            % error_measure(ii,jj,kk) = L2Norm;

            elapsedTime = toc;
            estimatedTotalTime = elapsedTime / prog;
            remainingTime = estimatedTotalTime - elapsedTime;
            
            if remainingTime < 120
                waitbar(prog, h, sprintf('Progress: %1.2f%%, Estimated time remaining: %1.1f seconds', 100*prog,  max(remainingTime, 0) )); % Update progress
            else
                waitbar(prog, h, sprintf('Progress: %1.2f%%, Estimated time remaining: %1.1f minutes', 100*prog,  max(remainingTime/60, 0) )); % Update progress
            end

        end
    end
end
close(h);
toc;

%%

figure;
x = linspace(P_A_vals(1), P_A_vals(end), length(P_A_vals));
y = linspace(P_D_vals(1), P_D_vals(end), length(P_D_vals));
[X, Y] = meshgrid(P_A_vals, P_D_vals);
cmap = myColour3Gradient(255,yellow_colour,red_colour,purple_colour);


% for kk=1:length(C_sens_vals)
for kk=1:length(C_sens_vals)

    subplot(1,3,kk)
    hold on;
    contour(X, Y, PDE_Activation(:,:,kk)', ':' ,'LineWidth',3)

    colormap(cmap);
    colorbar;

    contour(X, Y, ABM_Activation(:,:,kk)','-' ,'LineWidth',2)
    % surf(X, Y, ABM_Activation(:,:,kk)')
    title("$\chi = " + num2str(C_sens_vals(kk)) + "$",'FontSize',18,'Interpreter','latex');

    % legend_string = [legend_string; ["$\chi = " + num2str(C_sens_vals(kk)) + "$"]];

    xlabel('$\alpha_+$','FontSize',18,'Interpreter','latex');
    ylabel('$\alpha_-$','FontSize',18,'Interpreter','latex');
    
    xticks([0:0.5:1])
    xticklabels([0:0.5:1])
    yticks([0:0.5:1])
    yticklabels([0:0.5:1])

end

%%

figure;
x = linspace(P_A_vals(1), P_A_vals(end), length(P_A_vals));
y = linspace(P_D_vals(1), P_D_vals(end), length(P_D_vals));
[X, Y] = meshgrid(P_A_vals, P_D_vals);
cmap = myColour3Gradient(255,yellow_colour,red_colour,purple_colour);


% for kk=1:length(C_sens_vals)
for kk=1:length(C_sens_vals)

    subplot(1,2,1)
    % contour(X, Y, PDE_FullActivation(:,:,kk)'/ModelParams.activatedAge, ':' ,'LineWidth',3)
    surf(X, Y, PDE_FullActivation(:,:,kk)')
    view(2);

    colormap(cmap);
    colorbar;

    subplot(1,2,2)
    % contour(X, Y, ABM_FullActivation(:,:,kk)','-' ,'LineWidth',2)
    surf(X, Y, ABM_FullActivation(:,:,kk)')
    view(2);
    title("$\chi = " + num2str(C_sens_vals(kk)) + "$",'FontSize',18,'Interpreter','latex');
    xlabel('$\alpha_+$','FontSize',18,'Interpreter','latex');
ylabel('$\alpha_-$','FontSize',18,'Interpreter','latex');
    % legend_string = [legend_string; ["$\chi = " + num2str(C_sens_vals(kk)) + "$"]];
end



%%
% figure;
% hold on;
% 
% plot(P_A_vals,ABM_Activation,'-o','LineWidth',2)
% plot(P_A_vals,PDE_Activation,'-o','LineWidth',2)
% 
% hold off;

cmap = myColour3Gradient(length(C_sens_vals),yellow_colour,red_colour,purple_colour);

figure;
hold on;
x = linspace(P_A_vals(1), P_A_vals(end), length(P_A_vals));
y = linspace(P_D_vals(1), P_D_vals(end), length(P_D_vals));
[X, Y] = meshgrid(P_A_vals, P_D_vals);

legend_string = [];

% for kk=1:length(C_sens_vals)
for kk=1:length(C_sens_vals)
    % ModelParams.C_sens = C_sens_vals(kk);
    
    % activation_kk = (PDE_Activation(:,:,kk)-ABM_Activation(:,:,kk))';
    % h = contour(X, Y, activation_kk ,[-0.1,-0.1],'LineWidth',3,'color', cmap(kk,:));
    % surf(X, Y, activation_kk )
    
    error_measure_kk = error_measure(:,:,kk)';
    contour(X, Y, error_measure_kk ,[0.05, 0.05],':', 'LineWidth',3,'color', cmap(kk,:));
    contour(X, Y, error_measure_kk ,[0.10, 0.10],'--','LineWidth',3,'color', cmap(kk,:));
    contour(X, Y, error_measure_kk ,[0.15, 0.15],'-', 'LineWidth',3,'color', cmap(kk,:));

% h1 = contour(X, Y, error_measure_kk ,[1, 1], 'linestyle',':' ,'LineWidth',3,'color', cmap(kk,:));
    % h2 = contour(X, Y, error_measure_kk ,[0.25, 0.25], 'linestyle','--','LineWidth',3,'color', cmap(kk,:));
    % h3 = contour(X, Y, error_measure_kk ,[0.5, 0.5], 'linestyle','-','LineWidth',3,'color', cmap(kk,:));
    % 

    % surf(X, Y, error_measure_kk )
    colorbar

    legend_string = [legend_string; ["$\chi = " + num2str(C_sens_vals(kk)) + "$"]];
end
% surf(X, Y, (PDE_Activation-ABM_Activation)' )
xlabel('$\alpha_+$','FontSize',18,'Interpreter','latex');
ylabel('$\alpha_-$','FontSize',18,'Interpreter','latex');
% colorbar;
legend(legend_string,'Location','north','FontSize',16,'Interpreter','latex');
hold off;



%%

f = gcf;
f.Color = [1 1 1];
export_fig 2D_Comparision.png -m2.5



%%
function area_between = computeL1NormError(x1,y1,x2,y2)

    % Define two different x-meshes
    % x1 = 0:ModelParams.dt_PDE:ModelParams.T_final'; % First x-mesh
    % y1 = activatedCells_PDE'; % First curve
    % 
    % x2 = 0:ModelParams.dt_ABM:ModelParams.T_final'; % Second x-mesh
    % y2 = activatedCells_ABM'; % Second curve
    
    % Interpolate y2 onto x1's mesh
    y2_interp = interp1(x2, y2, x1, 'previous'); % Linear interpolation
    
    % Compute the absolute area between the curves
    % area_between = sqrt(trapz(x1, abs(y1 - y2_interp).^2));
    area_between = (trapz(x1, abs(y1 - y2_interp)));
    
    % disp(['Area between curves: ', num2str(area_between)]);

end

function L2Norm = computeL2NormError(x1,y1,x2,y2)

    % Define two different x-meshes
    % x1 = 0:ModelParams.dt_PDE:ModelParams.T_final'; % First x-mesh
    % y1 = activatedCells_PDE'; % First curve
    % 
    % x2 = 0:ModelParams.dt_ABM:ModelParams.T_final'; % Second x-mesh
    % y2 = activatedCells_ABM'; % Second curve
    
    % Interpolate y2 onto x1's mesh
    y2_interp = interp1(x2, y2, x1, 'previous'); % Linear interpolation
    
    % Compute the absolute area between the curves
    L2Norm = norm(y1-y2_interp);
    
    % disp(['Area between curves: ', num2str(area_between)]);

end