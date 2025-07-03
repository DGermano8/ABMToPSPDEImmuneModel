clear all;
close all;

%%
ModelParams = struct();
rndSeed = 106;

domainBoundary = struct();
domainBoundary.x_max =  60;
domainBoundary.y_max =  60;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 1000;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 0.5;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 16;
ModelParams.R_DC = 1.00;
ModelParams.numberOfClusters = 3;

% T-cells
ModelParams.P_A = 2 % Antigen gain rate
ModelParams.P_D = 0.5;  % Antigen loss rate
ModelParams.activatedAge  = 50;

% ModelParams.dx_ABM = 0.25;
% ModelParams.dy_ABM = 0.25;
% ModelParams.dt_ABM = 0.01;            % Time step
% ModelParams.NumbAgents = 5000;
% 
% ModelParams.dx_PDE = 0.25;
% ModelParams.dy_PDE = 0.25;
% ModelParams.dt_PDE = 0.01;

ModelParams.dx_ABM = 1/3;
ModelParams.dy_ABM = 1/3;
ModelParams.dt_ABM = 0.025;   % Time step
ModelParams.NumbAgents = 100000;

ModelParams.dx_PDE = 1/3;
ModelParams.dy_PDE = 1/3;
ModelParams.dt_PDE = 0.025;


if (ModelParams.p_move*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    ModelParams.C_sens*(4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) >= 1 || ...
    (4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2) > 1)
    clc;
    disp('Choose smaller dx and dt for ABM')
else
    clc;
    disp('Nice params!')
    disp([' ps = ', num2str((4*ModelParams.dt_ABM)/(ModelParams.dx_ABM.^2))] )

end

%%
tic;

numb_reps = 100;

poolobj = gcp("nocreate"); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(4);
end

iter=0;

activatedCells_PDE =   zeros(ModelParams.T_final/ModelParams.dt_PDE+1,numb_reps);
activatedCells_T_PDE = zeros(ModelParams.T_final/ModelParams.dt_PDE+1,numb_reps);
h = waitbar(0, 'Processing...'); % Create progress bar

parfor jj=1:numb_reps

    rndSeed = jj;
    [U, A, U_a, A_a, dc1_dx, dc2_dx, dc1_dy, dc2_dy, Ic, A_Ic, Inds_Ic, params_PDE] = PDE_SetUp(rndSeed,domainBoundary,ModelParams);

    
    tmp_activated = zeros(ModelParams.T_final/ModelParams.dt_PDE+1,1);
    tmp_total_activated = zeros(ModelParams.T_final/ModelParams.dt_PDE+1,1);
    for n = 0:params_PDE.nt
        
        % 
        % [walker_positions, walker_activation, U_abm, C_abm, Ap_abm, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
        % activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*length(walker_activation));
        % activatedCells_T_ABM(n+1) = nnz(walker_activation==params_ABM.activatedAge)/length(walker_activation);
        % 
        % 
        [U, A, U_a, A_a] = computePDEModel(U, A, U_a, A_a, Ic, Inds_Ic, A_Ic,  dc1_dx, dc2_dx, dc1_dy, dc2_dy, params_PDE);
        tmp = (A+A_a);
        tmp_activated(n+1) = ( sum( tmp(:) ) )/(ModelParams.activatedAge);
        tmp_total_activated(n+1) = ( sum( A_a(:) ) )/(ModelParams.activatedAge);
        
    end
    activatedCells_PDE(:,jj) = tmp_activated(:);
    activatedCells_T_PDE(:,jj) = tmp_total_activated(:);

    jj
    % iter = iter + 1;
    % prog = (iter)/(numb_reps);
    % 
    % elapsedTime = toc;
    % estimatedTotalTime = elapsedTime / prog;
    % remainingTime = estimatedTotalTime - elapsedTime;
    % 
    % if remainingTime < 120
    %     waitbar(prog, h, sprintf('Progress: %1.2f%%, Estimated time remaining: %1.1f seconds', 100*prog,  max(remainingTime, 0) )); % Update progress
    % else
    %     waitbar(prog, h, sprintf('Progress: %1.2f%%, Estimated time remaining: %1.1f minutes', 100*prog,  max(remainingTime/60, 0) )); % Update progress
    % end
end


%%

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];
white_colour = 1/255*[217 217 217];

plot_vet = 1:(100):(ModelParams.T_final/ModelParams.dt_PDE+1);
time_vect = ModelParams.dt_PDE*plot_vet;

figure;
subplot(2,2,1)
hold on;
for jj=1:numb_reps
    plot(time_vect, activatedCells_PDE(plot_vet,jj), Color=white_colour)
end
axis([0 ModelParams.T_final 0 1])
ylabel('Proportion of activation', 'FontSize',14, 'Interpreter','latex');
xlabel('Time', 'FontSize',14, 'Interpreter','latex');

subplot(2,2,3)
hold on;
for jj=1:numb_reps
    plot(time_vect, activatedCells_T_PDE(plot_vet,jj), Color=white_colour)
end
axis([0 ModelParams.T_final 0 1])
ylabel('Proportion of activated T cells', 'FontSize',14, 'Interpreter','latex');
xlabel('Time', 'FontSize',14, 'Interpreter','latex');

subplot(2,2,2)
violinplot(activatedCells_PDE(end,:))
axis([0.5 1.5 0 1])
xticklabels([3])
xlabel('Mean DC clusters', 'FontSize',14, 'Interpreter','latex');

subplot(2,2,4)
violinplot(activatedCells_T_PDE(end,:))
axis([0.5 1.5 0 1])
xticklabels([3])
xlabel('Mean DC clusters', 'FontSize',14, 'Interpreter','latex');

MakeDark();
