x_max = 10;
dx = 0.1;

dt = 0.01;
n_max = 50000;

nx = x_max/dx+1;

max_A = 50;
loss_rate = 0.5;
gain_rate = 0.9;

n_rand = 2000;

D_move = 0.5;
p_move = D_move * dx^2/(2*dt);

A_abm = zeros(n_max,n_rand);
R_abm = zeros(n_max,n_rand);
R_abm(1,:) = randi(nx,n_rand,1);


A_pde = zeros(n_max,nx);

A_abm(1,:) = 0;
A_pde(1,:) = max_A;



Ic = ones(nx,1);
boundary_split = 0.5*nx;
Ic(round(boundary_split):end) = 0;


for ii = 2:n_max
    %% do the ABM
    % update space 
    x_i = R_abm(ii-1,:);
    Rad = rand(n_rand,1);

    prob_left =  p_move/2;
    prob_right = p_move/2;
    prob_stay = (1-p_move);
    
    x_i( Rad < prob_left) = x_i( Rad < prob_left) - 1;
    x_i( Rad > prob_left & Rad < prob_left + prob_right) = x_i( Rad > prob_left &  Rad < prob_left + prob_right) + 1;
    x_i(x_i < 0) = 0;
    x_i(x_i > nx) = nx;
    
    R_abm(ii,:) = x_i;
    
    % update state
    a_i = A_abm(ii-1,:);
    near_boundary = (x_i < boundary_split);
    prob_up = near_boundary' .* gain_rate .* dt .* ( 1 - a_i'./max_A);
    prob_down = (1-near_boundary') .* loss_rate .* dt .*a_i' .* (1 - (a_i'==max_A) );
    prob_stay = 1 - prob_up - prob_down;

       % Draw a random event
    Rad = rand(n_rand,1);
    
    a_i( Rad < prob_up) = a_i( Rad < prob_up) + 1;
    a_i( Rad > prob_up & Rad < prob_up + prob_down) = a_i( Rad > prob_up & Rad < prob_up + prob_down) - 1;
    a_i(a_i > max_A) = max_A;
    a_i(a_i < 0) = 0;

    A_abm(ii,:) = a_i;
    
    if mod(ii,250) == 0
        a_x = zeros(nx,1);
        for jj=1:nx
            a_x(jj) = sum(a_i(x_i == jj)) / (max_A*n_rand);
        end

        clf;
        hold on
        subplot(2,2,1)
        histogram(x_i,1:5:nx)
        axis([1 nx 0 1.25*5*n_rand/(nx)])

        subplot(2,2,2)
        histogram(a_i,1:max_A)
        
        subplot(2,2,3)
        plot(a_x)
        

        subplot(2,2,4)
        plot(sum(A_abm(1:ii,:),2)/(n_rand*max_A) , 'LineWidth',2)
        axis([1 n_max 0 1])

        sgtitle(['time = ', num2str(dt*ii)])
        drawnow;
        hold off;
    end    
    % if ii < 0.25*n_max
    %     A_pde(ii) = A_pde(ii-1) - dt*loss_rate*A_pde(ii-1);
    % 
    %     for jj=1:n_rand
    %         if rand < dt*loss_rate*A_abm(ii-1,jj)
    %            A_abm(ii,jj) = A_abm(ii-1,jj) - 1;
    %         else
    %             A_abm(ii,jj) = A_abm(ii-1,jj);
    %         end
    %     end
    % else
    %     A_pde(ii) = A_pde(ii-1) + dt*gain_rate*(1-A_pde(ii-1)/max_A);
    % 
    %     for jj=1:n_rand
    %         if rand < dt*gain_rate*(1-A_abm(ii-1,jj)/max_A)
    %            A_abm(ii,jj) = A_abm(ii-1,jj) + 1;
    %         else
    %             A_abm(ii,jj) = A_abm(ii-1,jj);
    %         end
    %     end
    % end

end
