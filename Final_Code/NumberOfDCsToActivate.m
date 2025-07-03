clear all;
% close all;

addpath(genpath('UserFunctions'))


purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];
white_colour = 1/255*[217 217 217];

cmap = myColour2Gradient(255, red_colour, blue_colour);

plot_abm = true;
plot_pde = true;

%%

ModelParams = struct();
rndSeed = 106;

domainBoundary = struct();
domainBoundary.x_max =  70;
domainBoundary.y_max =  70;
domainBoundary.x_min =  0;
domainBoundary.y_min =  0;

ModelParams.T_final = 2000;
ModelParams.p_move = 0.5;

% ChemoTaxis params
ModelParams.C_sens = 1;             % Taxis sensitivity coefficient
ModelParams.C_scaling = 0.001;

% DCs
ModelParams.NumDCs = 16;
ModelParams.R_DC = 1.00;
% ModelParams.numberOfClusters = 8;
ModelParams.numberOfClusters = 1;

% T-cells
% ModelParams.P_A = 0.4; % Antigen gain rate
ModelParams.P_A = 0.8; % Antigen gain rate
ModelParams.P_D = 0.1;  % Antigen loss rate
ModelParams.activatedAge  = 20;

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
ModelParams.NumbAgents = 2000;

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
% h = waitbar(0, 'Processing...'); % Create progress bar

%%
numb_reps = 8;

NumberDCsToActivate = zeros(ModelParams.NumDCs,numb_reps);

poolobj = gcp("nocreate"); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(4);
end

tic;
parfor jj = 1:numb_reps
    rndSeed = jj; 
    [walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM] = Discrete_SetUp(rndSeed,domainBoundary,ModelParams);
    activatedCells_ABM = zeros(params_ABM.nt+1,1);
    activatedCells_T_ABM = zeros(params_ABM.nt+1,1);
    
    for n = 0:params_ABM.nt
        
        [walker_positions, walker_activation, U_abm, C_abm, Ap_abm, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params_ABM);
        activatedCells_ABM(n+1) = sum(walker_activation(:))/(params_ABM.activatedAge*length(walker_activation));
        activatedCells_T_ABM(n+1) = nnz(walker_activation==params_ABM.activatedAge)/length(walker_activation);
    
    end
    jj
    
    NumberDCsToActivate_tmp = zeros(ModelParams.NumDCs,1);
    for ii=1:size(DCLingerTime,1)
    
        if walker_activation(ii) == ModelParams.activatedAge
            ii_NumberDCsToActivate = nnz(DCLingerTime(ii,:));
            if ii_NumberDCsToActivate > 0
                NumberDCsToActivate_tmp(ii_NumberDCsToActivate) = NumberDCsToActivate_tmp(ii_NumberDCsToActivate) + 1;
            end
        end
    end
    NumberDCsToActivate(:,jj) = NumberDCsToActivate_tmp(:)/ModelParams.NumbAgents;
end
toc;
%%


figure;


% subplot(2,2,4)
% Plot the histogram using bar
if ModelParams.P_A < 0.5
    bar(1:ModelParams.NumDCs, sum(NumberDCsToActivate,2)/numb_reps, 'FaceColor', red_colour);
else
    bar(1:ModelParams.NumDCs, sum(NumberDCsToActivate,2)/numb_reps  , 'FaceColor', blue_colour);
end
% Add labels and title
xlabel('Number of DCs', 'FontSize',12, 'Interpreter','latex');
ylabel('Proportion of T cells activated', 'FontSize',12, 'Interpreter','latex');
axis([0 ModelParams.NumDCs+1 0 1.05*max(sum(NumberDCsToActivate,2)/numb_reps)])

MakeDark();

%%
