% close all;
clear all;
clc;

videoFile = 'simulation_video_4.avi';  % Change to .mp4 for MP4 format
v = VideoWriter(videoFile, 'MPEG-4'); % Use 'Uncompressed AVI' for lossless
v.FrameRate = 100;  % Adjust frame rate as needed
open(v);

purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];

x_max = 10;
dx = 0.2;

dt = 0.01;
max_time = 500;
n_max = max_time/dt;

(4*dt)/(dx.^2)

%%

plot_times = [0:25:max_time]/dt;
plot_t = length(plot_times);

nx = x_max/dx+1;

max_A = 100;
loss_rate = 0.5;
gain_rate = 0.5;

n_rand = 100000;

D_move = 0.5;
p_move = D_move * (2*dt)/dx^2;
chi_sense = 0.5;
p_chem = chi_sense * (2*dt)/dx^2;
if( p_move  > 1)
    disp('NO')
end

a_i = zeros(1,n_rand);
x_i = (nx-1)*ones(n_rand,1);

A_pde = zeros(nx,1);
U_pde = zeros(nx,1);
U_pde(end) = 1;


% right boundary
% Ic = zeros(nx,1);
% boundary_split = 0.5*nx;
% Ic(round(boundary_split):end) = 1;

% left boundary
Ic = ones(nx,1);
boundary_split = 0.1*nx;
Ic(round(boundary_split):end) = 0;

C =  1-(0:dx:x_max)'/x_max;

activation_abm = zeros(n_max,1);
activation_pde = zeros(n_max,1);
activation_T_abm = zeros(n_max,1);

%% define components for Phenotype model
Lx = x_max; 
La = max_A; 
T = max_time;

Dx = dx; 
Da = 1; 
Dt = dt;

Nt = T/Dt;
Na = La/Da;
Nx = x_max/Dx;

D = D_move;
chi = chi_sense; 
Da_diffusion = 0.01;
chia = 1; % 'a'-chemotaxis sensitivity
rho_plus = gain_rate;
rho_minus = loss_rate;

% Create grid
x = linspace(0, Lx, Nx); 
a = linspace(0, La, Na); 
[X, A] = meshgrid(x,a);
X=X'; A=A';
C_ph = X/x_max;

t = linspace(0, T, Nt + 1);

% Initialize u (example initial condition)
u = zeros(Nx, Na);
u(end,1) = 1;

activation_proportion = zeros(Nt + 1,1);

% Define f(u, a) (example function)
Ic_ph = zeros(Nx,Na);
Ic_ph(1:round(0.1*Nx),:) = 1;
% Ic(:,:) = 1;



%%
a_x = zeros(nx,1);

pde_line = 3;
abm_line = 2.2;
figure;

c_extended = padarray(C, [1], 'replicate'); % Extend boundaries by mirroring
% Extract shifted matrices (inner region only)
c_shift_xp = c_extended(3:end); % Shift down
c_shift_xm = c_extended(1:end-2); % Shift up
% For chemotactic gradients
dc1_dx = (c_shift_xp - C) / dx;
dc2_dx = (C - c_shift_xm) / dx;

for ii = 2:n_max

    % ABM
    [a_i, x_i, activation_abm, activation_T_abm] =computeABM_1D(ii, a_i, x_i, n_rand,p_move, Da, nx, C, p_chem, boundary_split, gain_rate, loss_rate, max_A, dt, activation_abm, activation_T_abm);
    
    % Phenotype PDE
    [u, activation_proportion] = computePheno_1D(ii, u, Nx, Na, Dx, Da, Dt, C_ph,  chi, D, chia, A, a, activation_proportion, La, rho_plus,rho_minus, Ic_ph );

    %% Plot both

    if ismember(ii,plot_times) && true
        a_x = zeros(nx,1);
        for jj=0:(nx-1)

            if any(x_i == jj)
                a_x(jj+1) = sum(a_i(x_i == jj))/(nnz(x_i == jj));
            end
        end
        
        clf;
        subplot(2,2,1)
        hold on
        [N, edges] = histcounts(x_i,-0.5:1:(nx-0.5));
        plot(dx*(1:nx), N/n_rand, 'LineWidth',2)
        plot(x ,sum(u,2), ':',LineWidth=4)
        title('position density')

        subplot(2,2,2)
        hold on;
        plot(dx*(1:nx), N'/n_rand.*a_x/max_A, 'LineWidth',2)
        avg_a = (u*a') /(La*sum(u(:)));
        plot(x ,avg_a/da, ':',LineWidth=4)
        title('activation density')
        axis([0 dx*nx 0 max([0.001; avg_a])*1.1 ])

        
        subplot(2,2,3)
        % histogram(a_i,0:max_A+1)
        hold on;
        [N_A, A_edges] = histcounts(a_i,-Da:Da:(max_A+Da));
        plot(A_edges(2:end-1), N_A(2:end)/(Da*n_rand), 'LineWidth',2)
        plot(a ,sum(u,1)/Da, ':',LineWidth=4)
        axis([0 max_A 0 max([N_A/n_rand, sum(u,1)/Da])*1.1])
        title('distribution of activation')
        

        subplot(2,2,4)
        hold on;
        plot(dt*(1:ii), activation_abm(1:ii), 'LineWidth',2)
        plot(dt*(1:ii), activation_proportion(1:ii), ':',LineWidth=4)
        axis([dt dt*n_max 0 1])
        title('proportion of activation')


        sgtitle(['time = ', num2str(dt*ii), ' \alpha_+ = ', num2str(gain_rate), ' \alpha_- = ', num2str(loss_rate), ])
        % sgtitle(['PhenMass = ', num2str(sum(u(:)))])

        % MakeDark();
        drawnow;
        hold off;

    end

end

%%
f = gcf;
f.Color = [1 1 1];
% export_fig 1D_NoTaxis_NoLoss.png -m2.5


%%
function [a_i, x_i, activation_abm, activation_T_abm] =computeABM_1D(ii, a_i, x_i, n_rand,p_move, Da, nx, C, p_chem, boundary_split, gain_rate, loss_rate, max_A, dt, activation_abm, activation_T_abm)
    % update space 
    % x_i = R_abm(ii-1,:);
    Rad = rand(n_rand,1);

    prob_left =  p_move/2;
    prob_right = p_move/2;
    % prob_stay = (1-p_move);
    
    x_i( Rad < prob_left) = x_i( Rad < prob_left) - 1;
    x_i( Rad > prob_left & Rad < prob_left + prob_right) = x_i( Rad > prob_left &  Rad < prob_left + prob_right) + 1;
    x_i(x_i < 0) = 0;
    x_i(x_i >= nx) = nx-1;


    ids = x_i+1;
    ids_m1 = ids - 1; ids_m1(ids_m1 < 1) = 1;
    ids_p1 = ids + 1; ids_p1(ids_p1 >= nx) = nx;
    C_here = C(ids);
    C_left = C(ids_m1);
    C_right = C(ids_p1);
    
    prob_left =  max(0, (C_left-C_here) .*(p_chem/2) );
    prob_right = max(0, (C_right-C_here).*(p_chem/2) );
    % prob_stay = (1-prob_left - prob_right);

    % Draw a random event
    Rad = rand(n_rand,1);

    x_i( Rad < prob_left) = x_i( Rad < prob_left) - 1;
    x_i( Rad > prob_left & Rad < prob_left + prob_right) = x_i( Rad > prob_left &  Rad < prob_left + prob_right) + 1;
    x_i(x_i < 0) = 0;
    x_i(x_i >= nx) = nx-1;
    
    % update state
    % right boundary
    % near_boundary = (x_i > boundary_split)';
    % left boundary
    near_boundary = (x_i <= boundary_split)';
    
    prob_up = near_boundary .* (dt*gain_rate/Da) .* ( 1 - a_i/max_A);
    prob_down = (1-near_boundary) .* (dt*loss_rate/Da) .*a_i/max_A .* (1 - (a_i==max_A) );
    prob_stay = 1 - prob_up - prob_down;

   % Draw a random event
    Rad = rand(n_rand,1)';
    
    a_i( Rad < prob_up) = a_i( Rad < prob_up) + Da;
    a_i( Rad > prob_up & Rad < prob_up + prob_down) = a_i( Rad > prob_up & Rad < prob_up + prob_down) - Da;
    a_i(a_i > max_A) = max_A;
    a_i(a_i < 0) = 0;

    activation_abm(ii) = sum(a_i,2)/(n_rand*max_A);
    activation_T_abm(ii) = nnz(a_i == max_A)/n_rand;
end

function [u, activation_proportion] = computePheno_1D(ii, u, nx,na, dx, da, dt, C_m,  chix, Dx, chia, A, a, activation_proportion, amax, rho_plus,rho_minus, Ic)
    
    dCdx = (C_m(2:end, :) - C_m(1:(end-1), :)) / dx;

    flux_west_diff = Dx * (u(2:end, :) - u(1:end-1, :)) / dx;
    u_west = 0.5 * (u(2:end, :) + u(1:end-1, :)); % Central difference for u at face
    flux_west_chem = -chix .* u_west .* dCdx; % Note the index for gradCx
    
    flux_west = (flux_west_diff - flux_west_chem) / dx;
    flux_west = cat(1, zeros(1,na), flux_west);
    
    flux_east_diff = Dx * (u(2:end, :) - u(1:end-1, :)) / dx;
    u_east = 0.5 * (u(2:end, :) + u(1:end-1, :)); % Central difference for u at face
    flux_east_chem = -chix * u_east .* dCdx;   % Note the index for gradCx
    
    flux_east = (flux_east_diff - flux_east_chem) / dx;
    flux_east = cat(1, flux_east, zeros(1,na));

    % --- Phenotypic Advection Term ---
    a_ind_1 = [2:na, na];
    a_ind_1 = a_ind_1(:);
    a_plus_half_a = (A(:,a_ind_1) + A(:,:)) / 2;
    f_k_plus_half = (rho_plus .* (1 - a_plus_half_a./amax) .* Ic(:, :) - rho_minus .* a_plus_half_a./amax .* (1-Ic(:, :)) ) ;
    u_upwind_plus =  u(:, a_ind_1) .*(f_k_plus_half <= 0) + u(:, :) .* (f_k_plus_half > 0) ;
    
    a_ind_2 = [1 1:(na-1)];
    a_ind_2 = a_ind_2(:);
    a_minus_half_a = (A(:,:) + A(:,a_ind_2) ) / 2;
    f_k_minus_half = (rho_plus .* (1 - a_minus_half_a./amax) .* Ic(:, :) - rho_minus .* a_minus_half_a./amax .* (1-Ic(:, :)) ) ;
    u_upwind_minus = u(:, :) .*( f_k_minus_half <= 0) + u(:, a_ind_2) .* (f_k_minus_half > 0);

    flux_a_plus = chia .* f_k_plus_half .* u_upwind_plus;
    flux_a_minus = chia .* f_k_minus_half .* u_upwind_minus;

    % Boundary in a (no flux)
    flux_a_minus(:,1) = 0;
    % Boundary in a (no flux)
    flux_a_plus(:,end) = 0;
    adv_phenotype = -(flux_a_plus - flux_a_minus) / da;

    u = u + dt*(flux_east - flux_west + adv_phenotype);

    activation_proportion(ii) = sum((u*a') / sum(u(:)))/amax;

end


