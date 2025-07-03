clear all;


max_A = 50;
loss_rate = 0.5;
gain_rate = 0.9;

dt = 0.01;

n_max = 1000;

n_rand = 100;

A_abm = zeros(n_max,n_rand);
A_pde = zeros(n_max,1);

A_abm(1,:) = max_A;
A_pde(1) = max_A;


for ii = 2:n_max
    % if ii < 0.5*n_max
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
    if ii < 0.5*n_max
        A_pde(ii) = A_pde(ii-1) - dt*loss_rate;

        for jj=1:n_rand
            if rand < dt*loss_rate
               A_abm(ii,jj) = A_abm(ii-1,jj) - 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        end
    else
        A_pde(ii) = A_pde(ii-1) + dt*gain_rate;
    
        for jj=1:n_rand
            if rand < dt*gain_rate
               A_abm(ii,jj) = A_abm(ii-1,jj) + 1;
            else
                A_abm(ii,jj) = A_abm(ii-1,jj);
            end
        end
    end
end

figure;
plot(A_pde)
hold on;
plot(mean(A_abm,2))

%%



max_A = 10;
loss_rate = 0.1;
gain_rate = 0.9;

dt = 0.01;

n_max = 10000;

n_rand = 1000;

A_abm = zeros(n_max,n_rand);
A_pde = zeros(n_max,1);

A_abm(1,:) = 0;
A_pde(1) = 0;


for ii = 2:n_max
    A_pde(ii) = A_pde(ii-1) - 0.5*dt*loss_rate*A_pde(ii-1) + 0.5*dt*gain_rate*(1-A_pde(ii-1)/max_A);

    for jj=1:n_rand
        if rand < 0.5
            if rand < dt*loss_rate*A_abm(ii-1,jj)
               A_abm(ii,jj) = A_abm(ii-1,jj) - 1;
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
plot(A_pde)
hold on;
plot(mean(A_abm,2))



