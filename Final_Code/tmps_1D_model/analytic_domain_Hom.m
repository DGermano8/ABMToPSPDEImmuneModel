

max_A = 50;
loss_rate = 0.5;
gain_rate = 0.9;

dt = 0.01;

n_max = 20000;

n_rand = 2000;

A_abm = zeros(n_max,n_rand);
A_pde = zeros(n_max,1);

A_abm(1,:) = max_A;
A_pde(1) = max_A;


for ii = 2:n_max
    if ii < 0.25*n_max
        A_pde(ii) = A_pde(ii-1) - dt*loss_rate*A_pde(ii-1)/max_A;

        for jj=1:n_rand
            if rand < dt*loss_rate*A_abm(ii-1,jj)/max_A
               A_abm(ii,jj) = A_abm(ii-1,jj) - 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        end
    else
        A_pde(ii) = A_pde(ii-1) + dt*gain_rate*(1-A_pde(ii-1)/max_A);
    
        for jj=1:n_rand
            if rand < dt*gain_rate*(1-A_abm(ii-1,jj)/max_A)
               A_abm(ii,jj) = A_abm(ii-1,jj) + 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        end
    end
end

figure;
plot(A_pde,LineWidth=4)
hold on;
plot(mean(A_abm,2),'--',LineWidth=4)

%%



max_A = 100;
loss_rate = 0.1
gain_rate = 0.5;

dt = 0.01;

n_max = 10000;

n_rand = 1000;

A_abm = zeros(n_max,n_rand);
A_pde = zeros(n_max,1);
A_T_pde = zeros(n_max,1);

A_abm(1,:) = 0;
A_pde(1) = 0;

trans_rate = 1.0

for ii = 2:n_max
    
    loss_term = 0.5*loss_rate/max_A*A_pde(ii-1);
    gain_term = 0.5*gain_rate/max_A*(max_A-(A_pde(ii-1)+A_T_pde(ii-1)));
    act_term = trans_rate*dt*gain_rate*(A_pde(ii-1)/max_A^2);

    A_pde(ii) = A_pde(ii-1) + dt*(gain_term - loss_term - act_term);
    A_T_pde(ii) = A_T_pde(ii) + trans_rate*dt*gain_rate*(A_pde(ii-1)/max_A);

    for jj=1:n_rand
        if rand < 0.5
            if rand < dt*loss_rate*A_abm(ii-1,jj)/max_A
               A_abm(ii,jj) = A_abm(ii-1,jj) - 1*(1 - A_abm(ii-1,jj)==max_A);
               % A_abm(ii,jj) = A_abm(ii-1,jj) - 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        else
            if rand < dt*gain_rate*(1-A_abm(ii-1,jj)/max_A)
               A_abm(ii,jj) = A_abm(ii-1,jj) + 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        end
    end

  
end

figure;
plot(A_pde+A_T_pde,LineWidth=4)
hold on;
plot(mean(A_abm,2),'--',LineWidth=4)

%%
% this is the version where we try to do it as a sum

max_A = 5;
loss_rate = 0.1
gain_rate = 0.5;

dt = 0.01;

n_max = 10000;

n_rand = 1000;

A_abm = zeros(n_max,n_rand);
A_pde = zeros(n_max,max_A+1);

A_abm(1,:) = 0;
A_pde(1,1) = 1;

trans_rate = 1.0

for ii = 2:n_max

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PDE
    
    A_pde_old = A_pde;
    for jj = 1:(max_A+1)

        if jj==1
            
            activation_term_jj =  gain_rate/max_A*( (jj+1)/max_A - A_pde_old(ii-1,jj) );
            decay_term_jjp1 =     loss_rate/max_A*A_pde_old(ii-1,jj+1);

            A_pde(ii,jj) = A_pde_old(ii-1,jj) + dt*(decay_term_jjp1 - activation_term_jj) ;
            
        elseif jj == (max_A+1)
            
            activation_term_jjm1 = gain_rate/max_A*( jj/max_A - sum(A_pde_old(ii-1,1:(jj-1)),2) );

            A_pde(ii,jj) = A_pde_old(ii-1,jj) + dt*(activation_term_jj) ;

        else

            activation_term_jjm1 = gain_rate/max_A*( jj/max_A - sum(A_pde_old(ii-1,1:(jj-1)),2) );
            activation_term_jj =   gain_rate/max_A*( (jj+1)/max_A - sum(A_pde_old(ii-1,1:(jj)),2) );
            decay_term_jj =        loss_rate/max_A*A_pde_old(ii-1,jj);
            decay_term_jjp1 =      loss_rate/max_A*A_pde_old(ii-1,jj+1);

            A_pde(ii,jj) = A_pde_old(ii-1,jj) + dt*(decay_term_jjp1 - decay_term_jj + activation_term_jjm1 - activation_term_jj) ;

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ABM
    for jj=1:n_rand
        if rand < 0.5
            if rand < dt*loss_rate*A_abm(ii-1,jj)/max_A
               A_abm(ii,jj) = A_abm(ii-1,jj) - 1*(1 - A_abm(ii-1,jj)==max_A);
               % A_abm(ii,jj) = A_abm(ii-1,jj) - 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        else
            if rand < dt*gain_rate*(1-A_abm(ii-1,jj)/max_A)
               A_abm(ii,jj) = A_abm(ii-1,jj) + 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end

figure;
plot(sum(A_pde,2),LineWidth=4)
hold on;
plot(mean(A_abm,2),'--',LineWidth=4)

%%
