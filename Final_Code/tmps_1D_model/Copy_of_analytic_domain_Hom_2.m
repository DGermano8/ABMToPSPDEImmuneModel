% close all;
clc;



x_max = 10;
dx = 0.5;

dt = 0.05;
max_time = 2000;
n_max = max_time/dt;

% plot_times = [25 50 100 200 300 400 500 750 1000]/dt;
plot_times = [0:20:max_time]/dt;
plot_t = length(plot_times);

nx = x_max/dx+1;

max_A = 10;
loss_rate = 0.5;
gain_rate = 0.5;


% if loss_rate == 0;
%     activation_power = 1.0
% 
% else
%     activation_power = 1.5
% end
activation_power = 2;
n_rand = 50000;

D_move = 0.25;
p_move = D_move * (2*dt)/dx^2;
chi_sense = 0.25;
p_chem = chi_sense * (2*dt)/dx^2;
if( p_move  > 1)
    disp('NO')
end

a_i = zeros(1,n_rand);
% x_i = randi(nx,n_rand,1);
x_i = (nx-1)*ones(n_rand,1);

A_pde = zeros(nx,1);
A_T_pde = zeros(nx,1);
U_pde = zeros(nx,1);
U_pde(end) = 1;

% right boundary
% Ic = zeros(nx,1);
% boundary_split = 0.5*nx;
% Ic(round(boundary_split):end) = 1;

% left boundary
Ic = ones(nx,1);
boundary_split = 0.5*nx;
Ic(round(boundary_split):end) = 0;

C =  1-(0:dx:x_max)'/x_max;

activation_abm = zeros(n_max,1);
activation_pde = zeros(n_max,1);

Full_activation_abm = zeros(n_max,1);
Full_activation_pde = zeros(n_max,1);


purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];


% cmap_abm = flipud( myColour3Gradient(n_max/plot_t+1,  1/255*[224 236 244], 1/255*[158 188 218], 1/255*[136 86 167]) );
% cmap_pde = flipud( myColour3Gradient(n_max/plot_t+1,  1/255*[255 237 160], 1/255*[254 178 76],  1/255*[240 59 32])  );
% 
cmap_abm = ( myColour2Gradient(n_max/plot_t+1, 1/255*[230 202 92], 1/255*[82 204 118]) );
cmap_pde = ( myColour2Gradient(n_max/plot_t+1, 1/255*[255 102 100], 1/255*[160 167 255])  );

% figure;
% hold on;
% plot(0:dx:x_max, C, '--','LineWidth',2.5 , 'Color', green_colour)
% plot(0:dx:x_max, Ic, '-','LineWidth',2.5 , 'Color', yellow_colour)
% hold off;
% legend('$C$','$1_{{A}}$', 'Interpreter','latex', 'FontSize',14)
% xlabel('$x$', 'Interpreter','latex', 'FontSize',14)
% title('Tissue topology', 'Interpreter','latex', 'FontSize',16)
% axis([0 x_max -0.0 1.0])
%%
a_x = zeros(nx,1);

pde_line = 3;
abm_line = 2.2;
figure;

iter = 1;
subplot(1,3,1)
hold on
[N, edges] = histcounts(x_i,-0.5:1:(nx-0.5));
plot3(0:dx:x_max, 0*ones(length(0:dx:x_max),1), N/n_rand, '-','LineWidth',abm_line , 'Color', blue_colour)
plot3(0:dx:x_max, 0*ones(length(0:dx:x_max),1), U_pde , '--', 'LineWidth',pde_line, 'Color', red_colour)
% plot(1:1:nx, U_pde , '--', 'LineWidth',2, 'Color', 'k')
axis([0 nx 0 5/(nx)])

subplot(1,3,2)
hold on
plot3(0:dx:x_max, 0*ones(length(0:dx:x_max),1), N'/n_rand.*a_x/max_A, 'LineWidth',abm_line, 'Color', blue_colour)
plot3(0:dx:x_max, 0*ones(length(0:dx:x_max),1), A_T_pde+A_pde, '--','LineWidth',pde_line,   'Color', red_colour)

axis([0 x_max 1 n_max/plot_t+1 0 0.03])


c_extended = padarray(C, [1], 'replicate'); % Extend boundaries by mirroring
% Extract shifted matrices (inner region only)
c_shift_xp = c_extended(3:end); % Shift down
c_shift_xm = c_extended(1:end-2); % Shift up
% For chemotactic gradients
dc1_dx = (c_shift_xp - C) / dx;
dc2_dx = (C - c_shift_xm) / dx;

for ii = 2:n_max
    
    %% do the ABM
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

    % prob_up = near_boundary .* gain_rate .* dt .* ( 1 - a_i/max_A);
    % prob_down = (1-near_boundary) .* loss_rate .* dt .*a_i/max_A .* (1 - (a_i==max_A) );
    prob_up = near_boundary .* gain_rate .* dt ; 
    prob_down = (1-near_boundary) .* loss_rate .* dt .* (1 - (a_i==max_A) );
    prob_stay = 1 - prob_up - prob_down;

   % Draw a random event
    Rad = rand(n_rand,1)';
    
    a_i( Rad < prob_up) = a_i( Rad < prob_up) + 1;
    a_i( Rad > prob_up & Rad < prob_up + prob_down) = a_i( Rad > prob_up & Rad < prob_up + prob_down) - 1;
    a_i(a_i > max_A) = max_A;
    a_i(a_i < 0) = 0;

    % A_abm(ii,:) = a_i;

    activation_abm(ii) = sum(a_i,2)/(n_rand*max_A);
    Full_activation_abm(ii) = sum(a_i==max_A)/(n_rand);


    %% update PDE
    U_old = U_pde;
    u_extended = padarray(U_old, [1], 'replicate'); % Extend boundaries by mirroring
    u_shift_xp = u_extended(3:end); % Shift down
    u_shift_xm = u_extended(1:end-2); % Shift up
    uxx = (u_shift_xp - 2 * U_old + u_shift_xm) / (dx)^2;

    % For du_dx and du_dy terms
    du1_dx = (U_old + u_shift_xp) / 2;
    du2_dx = (U_old + u_shift_xm) / 2;
    % Taxis terms
    u_taxis = (du1_dx .* dc1_dx - du2_dx .* dc2_dx) / dx;

    U_pde = U_old + dt * (D_move * (uxx) - chi_sense*(u_taxis));

    A_old = A_pde;
    A_T_old = A_T_pde;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a_extended = padarray(A_old , [1], 'replicate'); % Extend boundaries by mirroring
    a_shift_xp = a_extended(3:end); % Shift down
    a_shift_xm = a_extended(1:end-2); % Shift up
    axx = (a_shift_xp - 2 * (A_old ) + a_shift_xm) / (dx)^2;

     % For du_dx and du_dy terms
    da1_dx = (A_old  + a_shift_xp) / 2;
    da2_dx = (A_old  + a_shift_xm) / 2;
    % Taxis terms
    a_taxis = (da1_dx .* dc1_dx - da2_dx .* dc2_dx) / dx;
    
    delay_time = 0;
    if loss_rate > 0
        delay_time = delay_time + max_A/(loss_rate);
    end
    if gain_rate > 0
        delay_time = delay_time + max_A/(gain_rate);
    end
    

    final_state_change = gain_rate/(max_A^activation_power)*Ic.*(A_old) * (ii*dt > delay_time);
    
    % final_state_change = gain_rate/(max_A^(1+1/(max_A)))*Ic.*(A_old);
    % final_state_change = gain_rate/(max_A)*Ic.*max(0,(A_old)./(U_old+A_old)); 

    A_pde = A_old + dt * ( D_move * (axx) - chi_sense*a_taxis + gain_rate/max_A*Ic.*(U_old - (A_old+A_T_old))  - loss_rate/max_A.*(1-Ic).*A_old - final_state_change );%

    a_extended = padarray(A_T_old , [1], 'replicate'); % Extend boundaries by mirroring
    a_shift_xp = a_extended(3:end); % Shift down
    a_shift_xm = a_extended(1:end-2); % Shift up
    axx = (a_shift_xp - 2 * (A_T_old ) + a_shift_xm) / (dx)^2;

    % For du_dx and du_dy terms
    da1_dx = (A_T_old  + a_shift_xp) / 2;
    da2_dx = (A_T_old  + a_shift_xm) / 2;
    % Taxis terms
    a_taxis = (da1_dx .* dc1_dx - da2_dx .* dc2_dx) / dx;

    A_T_pde = A_T_old  + dt * ( D_move * (axx) - chi_sense*a_taxis + final_state_change ) ;%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a_extended = padarray(A_old + A_T_old, [1], 'replicate'); % Extend boundaries by mirroring
    % a_shift_xp = a_extended(3:end); % Shift down
    % a_shift_xm = a_extended(1:end-2); % Shift up
    % axx = (a_shift_xp - 2 * (A_old + A_T_old) + a_shift_xm) / (dx)^2;
    % 
    %  % For du_dx and du_dy terms
    % da1_dx = (A_old + A_T_old + a_shift_xp) / 2;
    % da2_dx = (A_old + A_T_old + a_shift_xm) / 2;
    % % Taxis terms
    % a_taxis = (da1_dx .* dc1_dx - da2_dx .* dc2_dx) / dx;
    % 
    % delay_time = 0;
    % if loss_rate > 0
    %     delay_time = delay_time + max_A/loss_rate;
    % end
    % if gain_rate > 0
    %     delay_time = delay_time + max_A/gain_rate;
    % end
    % final_state_change = gain_rate/(max_A^1.5)*Ic.*(A_old) * (ii*dt > delay_time);
    % 
    % % final_state_change = gain_rate/(max_A^1.5)*Ic.*(A_old); 
    % 
    % A_pde = A_old + dt * ( D_move * (axx) - chi_sense*a_taxis + gain_rate/max_A*Ic.*(U_old - (A_old+A_T_old))  - loss_rate.*(1-Ic).*A_old - final_state_change );%
    % 
    % a_extended = padarray(A_T_old + A_old, [1], 'replicate'); % Extend boundaries by mirroring
    % a_shift_xp = a_extended(3:end); % Shift down
    % a_shift_xm = a_extended(1:end-2); % Shift up
    % axx = (a_shift_xp - 2 * (A_T_old + A_old) + a_shift_xm) / (dx)^2;
    % 
    % % For du_dx and du_dy terms
    % da1_dx = (A_T_old + A_old + a_shift_xp) / 2;
    % da2_dx = (A_T_old + A_old + a_shift_xm) / 2;
    % % Taxis terms
    % a_taxis = (da1_dx .* dc1_dx - da2_dx .* dc2_dx) / dx;
    % 
    % A_T_pde = A_T_old  + dt * ( D_move * (axx) - chi_sense*a_taxis + final_state_change ) ;%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    activation_pde(ii) = sum(A_T_pde+A_pde);
    Full_activation_pde(ii) = sum(A_T_pde);
    
    %% Plot both
    % if mod(ii,plot_t) == 0 ; %10/dt) == 0
    if ismember(ii,plot_times) %10/dt) == 0
        % pause(0.1)
        a_x = zeros(nx,1);
        for jj=0:(nx-1)
            % a_x(jj) = sum(a_i(x_i == jj)) / (max_A * nnz(x_i == jj));
            % a_x(jj) = sum(a_i(x_i == jj)) / (max_A);
            if any(x_i == jj)
                a_x(jj+1) = sum(a_i(x_i == jj))/(nnz(x_i == jj));
            end
        end
        if true
            clf;
            subplot(2,3,1)
            hold on
            [N, edges] = histcounts(x_i,-0.5:1:(nx-0.5));
            plot(1:nx,N/n_rand, 'LineWidth',2)
            plot(1:1:nx, U_pde , '--', 'LineWidth',2)
            % axis([0 nx 0 0.05])
    
            
            
            subplot(2,3,4)
            hold on;
            plot(N'/n_rand.*a_x/max_A, 'LineWidth',2)
            plot(A_T_pde+A_pde, '--','LineWidth',2)
    
            % plot(Ic, '--')
            % axis([0 nx 0 0.05])
    
            subplot(2,3,2)
            hold on;
            plot(dt*(1:n_max), activation_abm-Full_activation_abm, 'LineWidth',2)
            plot(dt*(1:n_max), activation_pde-Full_activation_pde, '--','LineWidth',2)
            
            subplot(2,3,5)
            hold on;
            plot(dt*(1:n_max), Full_activation_abm, 'LineWidth',2)
            plot(dt*(1:n_max), Full_activation_pde, '--','LineWidth',2)
            
            % axis([dt dt*n_max 0 1])


            subplot(2,3,3)
            histogram(a_i,1:max_A)


            sgtitle(['time = ', num2str(dt*ii)])
            drawnow;
            hold off;

        end

        if false
            iter = iter + 1;
            subplot(1,3,1)
            hold on
            [N, edges] = histcounts(x_i,-0.5:1:(nx-0.5));
            plot3(0:dx:x_max, ii*dt*ones(length(0:dx:x_max),1), N/n_rand, '-','LineWidth',abm_line , 'Color', blue_colour)
            plot3(0:dx:x_max, ii*dt*ones(length(0:dx:x_max),1), U_pde , '--', 'LineWidth',pde_line, 'Color', red_colour)
            % plot(1:1:nx, U_pde , '--', 'LineWidth',2, 'Color', 'k')
            axis([0 x_max 0 max_time 0 0.05])
            
            view(60,20)
            grid on
            % title('Cell density')
            xlabel('$x$', 'Interpreter','latex', 'FontSize',14)
            ylabel('$t$', 'Interpreter','latex', 'FontSize',14)
            zlabel('$u$', 'Interpreter','latex', 'FontSize',14)
            title('Local T cell density', 'Interpreter','latex', 'FontSize',16)

            subplot(1,3,2)
            hold on
            plot3(0:dx:x_max, ii*dt*ones(length(0:dx:x_max),1), N'/n_rand.*a_x/max_A, 'LineWidth',abm_line, 'Color', blue_colour)
            plot3(0:dx:x_max, ii*dt*ones(length(0:dx:x_max),1), A_T_pde+A_pde, '--','LineWidth',pde_line,   'Color', red_colour)
            
            axis([0 x_max 0 max_time 0 0.03])
            
            view(60,20)
            grid on
            % title('Antigen density')
            xlabel('$x$', 'Interpreter','latex', 'FontSize',14)
            ylabel('$t$', 'Interpreter','latex', 'FontSize',14)
            zlabel('$a + \hat{a}$', 'Interpreter','latex', 'FontSize',14)
            title('Local antigen density', 'Interpreter','latex', 'FontSize',16)

            drawnow;
        end

    end



end

 %%
% subplot(2,3,3)
% hold on;
% plot(dt*(1:n_max), activation_abm-Full_activation_abm,'-', 'LineWidth', 3 , 'Color', blue_colour)
% plot(dt*(1:n_max), activation_pde-Full_activation_pde, ':','LineWidth',pde_line,  'Color', red_colour)
% 
% title('Proportion of activated T cells', 'Interpreter','latex', 'FontSize',16)
% ylabel('T cell activation', 'Interpreter','latex', 'FontSize',14)
% xlabel('$t$', 'Interpreter','latex', 'FontSize',14)
% axis([dt dt*n_max 0 1])
% 
% subplot(2,3,6)
% hold on;
% plot(dt*(1:n_max), Full_activation_abm, 'LineWidth', 3 , 'Color', blue_colour)
% plot(dt*(1:n_max), Full_activation_pde, ':','LineWidth',pde_line,  'Color', red_colour)
% 
% 
% subplot(1,3,1)
% xticks([0 5 10])
% xticklabels([10 5 0])
% 
% subplot(1,3,2)
% xticks([0 5 10])
% xticklabels([10 5 0])

%%
f = gcf;
f.Color = [1 1 1];
% export_fig 1D_NoTaxis_NoLoss.png -m2.5


