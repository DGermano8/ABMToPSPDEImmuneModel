close all;
clc;
figure;

x_max = 10;
dx = 0.2;
dy = 0.2;

dt = 0.01;
max_time = 500;
n_max = max_time/dt;

nx = x_max/dx+1;

max_A = 50;
loss_rate = 0.0;
gain_rate = 0.5;

n_rand = 2000;

D_move = 0.05;
p_move = D_move * (2*dt)/dx^2;
if( p_move  > 1)
    disp('NO')
end

a_i = zeros(1,n_rand);
% x_i = randi(nx,n_rand,1);
x_i = 0*ones(n_rand,1);



A_pde = zeros(nx,1);
A_T_pde = zeros(nx,max_A);
U_pde = zeros(nx,1);
U_pde(1) = 1;

% right boundary
% Ic = zeros(nx,1);
% boundary_split = 0.5*nx;
% Ic(round(boundary_split):end) = 1;

% left boundary
Ic = ones(nx,1);
boundary_split = 2*nx;
Ic(round(boundary_split):end) = 0;

activation = zeros(n_max,1);
for ii = 2:n_max
    % do the ABM
    % update space 
    % x_i = R_abm(ii-1,:);
    Rad = rand(n_rand,1);

    prob_left =  p_move/2;
    prob_right = p_move/2;
    prob_stay = (1-p_move);
    
    x_i( Rad < prob_left) = x_i( Rad < prob_left) - 1;
    x_i( Rad > prob_left & Rad < prob_left + prob_right) = x_i( Rad > prob_left &  Rad < prob_left + prob_right) + 1;
    x_i(x_i < 0) = 0;
    x_i(x_i > nx) = nx;
    
    % R_abm(ii,:) = x_i;
    
    % update state
    % a_i = A_abm(ii-1,:);
    % right boundary
    % near_boundar*y = (x_i > boundary_split)';
    % left boundary
    near_boundary = (x_i <= boundary_split)';

    prob_up = near_boundary .* gain_rate .* dt .* ( 1 - a_i./max_A);
    prob_down = (1-near_boundary) .* loss_rate .* dt .*a_i .* (1 - (a_i==max_A) );
    prob_stay = 1 - prob_up - prob_down;

   % Draw a random event
    Rad = rand(n_rand,1)';
    
    a_i( Rad < prob_up) = a_i( Rad < prob_up) + 1;
    a_i( Rad > prob_up & Rad < prob_up + prob_down) = a_i( Rad > prob_up & Rad < prob_up + prob_down) - 1;
    a_i(a_i > max_A) = max_A;
    a_i(a_i < 0) = 0;

    % A_abm(ii,:) = a_i;

    activation(ii) = sum(a_i,2)/(n_rand*max_A);



    % update PDE
    U_old = U_pde;
    u_extended = padarray(U_old, [1], 'replicate'); % Extend boundaries by mirroring
    u_shift_xp = u_extended(3:end); % Shift down
    u_shift_xm = u_extended(1:end-2); % Shift up
    uxx = (u_shift_xp - 2 * U_old + u_shift_xm) / (dx)^2;

    U_pde = U_old + dt * (D_move * (uxx));

    A_old = A_pde;
    a_extended = padarray(A_old, [1], 'replicate'); % Extend boundaries by mirroring
    a_shift_xp = a_extended(3:end); % Shift down
    a_shift_xm = a_extended(1:end-2); % Shift up
    axx = (a_shift_xp - 2 * A_old + a_shift_xm) / (dx)^2;
    
    A_pde = A_old + dt * ( D_move * (axx) + gain_rate*Ic.*(1 - A_old./max_A ).*U_old/dx );%
       
    if mod(ii,10/dt) == 0
        a_x = zeros(nx,1);
        for jj=1:nx
            % a_x(jj) = sum(a_i(x_i == jj)) / (max_A * nnz(x_i == jj));
            % a_x(jj) = sum(a_i(x_i == jj)) / (max_A);
            if any(x_i == jj)
                a_x(jj) = sum(a_i(x_i == jj))/(nnz(x_i == jj));
            end
        end

        clf;
        subplot(2,2,1)
        hold on
        [N, edges] = histcounts(x_i,-0.5:1:(nx+0.5));
        plot(0:1:nx,N/n_rand, 'LineWidth',2)
        plot(1:1:nx, U_pde, 'LineWidth',2)
        axis([0 nx 0 5/(nx)])

        subplot(2,2,2)
        histogram(a_i,1:max_A)
        
        subplot(2,2,3)
        hold on;
        plot(a_x/max_A, 'LineWidth',2)
        plot(A_pde/max_A, 'LineWidth',2)

        plot(Ic, '--')
        % axis([0 nx 0 1.1])

        subplot(2,2,4)
        plot(dt*(1:n_max), activation, 'LineWidth',2)
        axis([dt dt*n_max 0 1])

        sgtitle(['time = ', num2str(dt*ii)])
        drawnow;
        hold off;
    end



end
