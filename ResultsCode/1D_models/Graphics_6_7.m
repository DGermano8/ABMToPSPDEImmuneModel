% close all;
% clear all;
% clc;

fps = 10;
fileName = 'test_1';



purple_colour = 1/255*[131, 96, 150];
red_colour = 1/255*[237, 123, 123];
yellow_colour = 1/255*[240, 184, 110];
blue_colour = 1/255*[106 123 237];
green_colour = 1/255*[77 150 138];
back_colour = 1/255*[56 56 56];
font_colour = 1/255*[217 217 217];
% rng(2);

%%
clc;

% Initialise PDE paramaters:

dt = 0.01;
max_time = 600;

n_max = max_time/dt + 1;
plot_times = [0:10:max_time]/dt;

cmap = myColour3Gradient(255, yellow_colour, red_colour, purple_colour );
% cmap = myColour2Gradient(length(plot_times), yellow_colour, purple_colour );
c_it = 0;

x_max = 10;
dx = 0.2;
nx = x_max/dx+1;

(2*dt)/(dx.^2)

da = 1; 
max_A = 100;
na = max_A/da + 1;


D_move = 0.5;
chi_sense = 0.5;
chia = 1; % 'a'-chemotaxis sensitivity

loss_rate = 1;
gain_rate = 1;

p_move = D_move * (2*dt)/dx^2;
p_chem = chi_sense * (2*dt)/dx^2;
p_gain = (dt*gain_rate/da);
p_loss = (dt*loss_rate/da);
n_rand = 10000;

%%%%%

impose_boundary = 0.64;
impose_boundary = 0.4;



% Initialise ABM
Ic = zeros(nx,1);
boundary_split = impose_boundary*nx;
Ic(round(boundary_split):end) = 1;
C = (0:dx:x_max)'/x_max;

activation_abm = zeros(n_max,1);

a_i = zeros(1,n_rand);
x_i = zeros(n_rand,1);

% Initialise PDE
Ic_ph = zeros(nx,na);
Ic_ph(round(impose_boundary*nx):end,:) = 1;
x = linspace(0, x_max, nx); 
a = linspace(0, max_A, na); 
[X, A] = meshgrid(x,a);
X=X'; A=A';
C_ph = 1-X/x_max;

% figure;
% hold on;
% plot(x,1-C_ph(:,1), '--', 'linewidth', 2,'Color', green_colour, 'DisplayName', '$C$')
% plot(x,Ic_ph(:,1), '-', 'linewidth', 2,'Color', yellow_colour, 'DisplayName', '$1_A$')
% legend('FontSize',12 ,'Interpreter','latex','Color',back_colour);
% hold off;
% f1 = gcf;
% f1.Position = [271 348 386 275];
% MakeDark();
% 
% export_fig 1D_SingleDC.png -m2.5;

activation_proportion = zeros(n_max,1);

u = zeros(nx, na);
u(1,1) = 1;


abm_line = 2.5;
pde_line = 2.5;





a_an = a;
C_an = C;
A_Ic_an = Ic_ph(:,1);

u_SS_sol = exp(chi_sense/D_move .* C_an);
u_SS_sol = u_SS_sol/sum(u_SS_sol(:));

a_SS_sol = A_Ic_an .* u_SS_sol;
Prob_A = sum(a_SS_sol(:))
Prob_Ac = 1 - Prob_A;

% sigma_val = 0.125;
sigma_val = 3 * (D_move) / (impose_boundary*(1-impose_boundary)*x_max*x_max);

alpha_val = sigma_val *  Prob_A  * max_A / loss_rate; 
beta_val =  sigma_val *  Prob_Ac * max_A / gain_rate;
U_scale = (1./beta(alpha_val,beta_val)) * (1/(max_A).^(alpha_val+beta_val)) .* (alpha_val+beta_val)./(alpha_val+ (gain_rate/loss_rate) * beta_val);
U_A =                                     a_an .*                           a_an.^(alpha_val-1) .*(max_A - a_an).^(beta_val-1);
U_Ac =  (gain_rate/loss_rate) .* (max_A - a_an) .* a_an.^(alpha_val-1) .*(max_A - a_an).^(beta_val-1);
kappa_1 = alpha_val
kappa_2 = beta_val

mean_act = (max_A * kappa_1)/(kappa_1+kappa_2+1) * ((1/(kappa_1 + (gain_rate/loss_rate)*kappa_2)) + 1)
mean_act_level = mean_act/max_A
u_of_a = U_scale.*(U_A+U_Ac)*da;



%%

CellDensity_abm_hist = zeros(length(plot_times), length(-0.5:1:(nx-0.5))-1);
AntigenDensity_abm_hist = zeros(length(plot_times), length(-0.5:1:(nx-0.5))-1);
AntigenDist_abm_hist = zeros(length(plot_times),length(-da:da:(max_A+da))-2);

CellDensity_pde_hist = zeros(length(plot_times), length(x));
AntigenDensity_pde_hist = zeros(length(plot_times), length(x));
AntigenDist_pde_hist = zeros(length(plot_times),length(a));


%%

f = figure;
f.Units = 'pixels';
f.Position = [106 360 1200 400];
pause(0.1);


vidfile = VideoWriter(fileName,'MPEG-4');
vidfile.FrameRate = fps;
open(vidfile);
numberOfSteps = 0;



for ii = 0:n_max

    % ABM
    [a_i, x_i, activation_abm] =computeABM_1D(ii, a_i, x_i, n_rand,p_move, da, nx, C, p_chem, boundary_split, p_gain, p_loss, max_A, dt, activation_abm);
    
    % Phenotype PDE
    [u, activation_proportion] = computePheno_1D(ii, u, nx, na, dx, da, dt, C_ph,  chi_sense, D_move, chia, A, a, activation_proportion, max_A, gain_rate,loss_rate, Ic_ph );


    % Plot both
    if ismember(ii,plot_times) && true
        pause(0.01);
        clf;
        c_it = c_it + 1;

        time_ii = dt*ii;
        
        a_x = zeros(nx,1);
        for jj=0:(nx-1)

            if any(x_i == jj)
                a_x(jj+1) = sum(a_i(x_i == jj))/(nnz(x_i == jj));
            end
        end
        
        % store history of params:
        [N, edges] = histcounts(x_i,-0.5:1:(nx-0.5));
        [N_A, A_edges] = histcounts(a_i,-0.5*da:da:(max_A+0.5*da));
        avg_a = (u*a') /(da*max_A*sum(u(:)));


        CellDensity_abm_hist(c_it,:) = N/n_rand;
        AntigenDensity_abm_hist(c_it,:) = N'/n_rand.*a_x/(da*max_A);
        AntigenDist_abm_hist(c_it,:) = N_A(1:end)/(da*n_rand);
        
        CellDensity_pde_hist(c_it,:) = sum(u,2);
        AntigenDensity_pde_hist(c_it,:) = avg_a;
        AntigenDist_pde_hist(c_it,:) = sum(u,1)/da;



        subplot(1,2,1)
        hold on

        plot( dx*(1:nx), N/n_rand,'-','LineWidth',abm_line, 'Color',purple_colour)
        plot(x ,sum(u,2),       '--','LineWidth',pde_line, 'Color',yellow_colour)
        plot(x,u_SS_sol,':','LineWidth',2.5, 'color',blue_colour)

        verts = [x_max*impose_boundary 0; x_max 0; x_max 1; x_max*impose_boundary 1];
        patch(verts(:,1),verts(:,2),font_colour,'Facealpha',0.1,'EdgeColor','none')
        plot([impose_boundary*x_max impose_boundary*x_max], [-1 1], ':', 'LineWidth',1.0, 'Color', font_colour)
        
        lg_str = [{'ABM', 'PS-PDE', 'Approx.'}];
        lg = legend(lg_str, 'FontSize',12 ,'Interpreter','latex','location','southeast', 'Color',back_colour);
        title('T cell density', 'FontSize',14 ,'Interpreter','latex')
        axis([0 x_max 0 0.04])
        xlabel('x', 'FontSize',12 ,'Interpreter','latex')
        yticks([0:0.01:0.04])
        
        MakeDark();


        subplot(1,2,2)
        hold on;
        plot(A_edges(1:end-1)+0.5*da, N_A(1:end)/(n_rand), '-', 'LineWidth',abm_line, 'Color',purple_colour)
        plot( a ,sum(u,1),                                   '--' ,'LineWidth',pde_line, 'Color',yellow_colour)
        plot(a_an,u_of_a,':','LineWidth',2.5, 'color',blue_colour)

        axis([0 max_A 0 0.05])
        title('Antigen distribution', 'FontSize',14 ,'Interpreter','latex')
        xlabel('a', 'FontSize',12 ,'Interpreter','latex')
        hold off;
        yticks([0:0.01:0.05])

        MakeDark();


        % subplot(1,3,3)
        % hold on;
        % % scatter(-1,-1,'Marker','none')
        % % scatter(-1,-1,'Marker','none')
        % % scatter(-1,-1,'Marker','none')
        % plot(dt*(1:(ii)), activation_abm(1:ii),        '-' ,'LineWidth',abm_line, 'Color',purple_colour)
        % plot(dt*(1:(ii)), activation_proportion(1:(ii)), '--','LineWidth',pde_line, 'Color',yellow_colour)
        % % plot([0 max_time], [mean_act_level mean_act_level],':','LineWidth',2.5, 'color',blue_colour)
        % 
        % axis([dt dt*n_max 0 1])
        % title('Proportion of T cell activation', 'FontSize',14 ,'Interpreter','latex')
        % xlabel('t', 'FontSize',12 ,'Interpreter','latex')
        % hold off;

        MakeDark();
        sgtitle(['Time = ' num2str(time_ii)], 'FontSize',14 ,'Interpreter','latex')
        drawnow;

        numberOfSteps = numberOfSteps + 1;
        F(numberOfSteps) = getframe(f);
        writeVideo(vidfile,F(numberOfSteps));

    end

end

pause(0.5);
close(vidfile)


% lg_str = [{['$A_{max}=' , num2str(max_A),'$'] , ['$\kappa_{1} \approx' , num2str(round(kappa_1,2)) ,'$'], ...
%            ['$ \kappa_{2} \approx' ,num2str(round(kappa_2,2)),'$'], ...
%            'ABM', 'PS-PDE','Approx.'}];
% 
% 
% lg = legend(lg_str, 'FontSize',16 ,'Interpreter','latex',NumColumns=2);
% % lg.Position(1:2) = [.7 .2];
% 
% f = gcf;
% f.Position = [1 469 1512 295];
% f.Color = [1 1 1];
% lg.Position(1:2) = [.68 .2];

%%















%%
function [a_i, x_i, activation_abm] =computeABM_1D(ii, a_i, x_i, n_rand,p_move, Da, nx, C, p_chem, boundary_split, p_gain, p_loss, max_A, dt, activation_abm)
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
    near_boundary = (x_i > boundary_split)';
    % left boundary
    % near_boundary = (x_i <= boundary_split)';
    
    prob_up = near_boundary .* (p_gain) .* ( 1 - (a_i+Da)/max_A);
    prob_down = (1-near_boundary) .* (p_loss) .*(a_i-Da)/max_A .* (1 - (a_i==max_A) );
    prob_stay = 1 - prob_up - prob_down;

   % Draw a random event
    Rad = rand(n_rand,1)';
    
    a_i( Rad < prob_up) = a_i( Rad < prob_up) + Da;
    a_i( Rad > prob_up & Rad < prob_up + prob_down) = a_i( Rad > prob_up & Rad < prob_up + prob_down) - Da;
    a_i(a_i > max_A) = max_A;
    a_i(a_i < 0) = 0;

    activation_abm((ii+1)) = sum(a_i,2)/(n_rand*max_A);
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

    activation_proportion((ii+1)) = sum((u*a') / sum(u(:)))/amax;

end


