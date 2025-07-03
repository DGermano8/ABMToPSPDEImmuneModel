
function [walker_positions, walker_activation, U, C, A, DCLingerTime] = computeABMModel_Vectorized(walker_positions, walker_activation, C, DCLingerTime, DC_model, params)
    Ic = DC_model.Ic;
    % DC_id = DC_model.DC_id;
    BoundaryDC = DC_model.BoundaryDC;

    U = zeros(params.Ny, params.Nx);
    A = zeros(params.Ny, params.Nx);
    
        % Get the walker's current position
        wx = walker_positions(:, 1);
        wy = walker_positions(:, 2);

        %% based on Random Motion
        id_x = round(wx/params.dx + 1);
        id_y = round(wy/params.dy + 1);

        % Compute next steps
        id_x_p1 = id_x+1;
        id_x_m1 = id_x-1;
        id_y_p1 = id_y+1;
        id_y_m1 = id_y-1;
        id_x_p1(id_x_p1 > params.Nx) =  params.Nx;
        id_x_m1(id_x_m1 < 1) =  1;
        id_y_p1(id_y_p1 > params.Ny) =  params.Ny;
        id_y_m1(id_y_m1 < 1) =  1;

        left_indices = sub2ind(size(Ic), id_y, id_x_m1);
        right_indices = sub2ind(size(Ic), id_y, id_x_p1);
        down_indices = sub2ind(size(Ic), id_y_m1, id_x);
        up_indices = sub2ind(size(Ic), id_y_p1, id_x);

        % Compute movement probabilities
        prob_left =  max(0, (1-Ic(left_indices))*(params.p_move/4) );
        prob_right = max(0, (1-Ic(right_indices))*(params.p_move/4) );
        prob_down =  max(0, (1-Ic(down_indices))*(params.p_move/4) );
        prob_up =    max(0, (1-Ic(up_indices))*(params.p_move/4) );
        prob_stay = (1-params.p_move)*ones(params.num_walkers,1);

        % Normalize probabilities (total probability of moving = 1)
        total_prob = prob_left + prob_right + prob_down + prob_up + prob_stay;
        prob_left = prob_left ./ total_prob;
        prob_right = prob_right ./ total_prob;
        prob_down = prob_down ./ total_prob;
        prob_up = prob_up ./ total_prob;
        prob_stay = prob_stay./total_prob;

        p1 = prob_left;
        p2 = prob_left + prob_right;
        p3 = prob_left + prob_right + prob_down;
        p4 = prob_left + prob_right + prob_down + prob_up;
        

        Rad=rand(params.num_walkers,1);

        % id_x(Rad < p1 & Rad > 0)  = id_x(Rad < p1 & Rad > 0)  - 1;
        % id_x(Rad < p2 & Rad > p1) = id_x(Rad < p2 & Rad > p1) + 1;
        % id_y(Rad < p3 & Rad > p2) = id_y(Rad < p3 & Rad > p2) - 1;
        % id_y(Rad < p4 & Rad > p3) = id_y(Rad < p4 & Rad > p3) + 1;
        % 
        % id_x(id_x > params.Nx) =  params.Nx;
        % id_x(id_x < 1) =  1;
        % id_y(id_y > params.Ny) =  params.Ny;
        % id_y(id_y < 1) =  1;

        id_x(Rad < p1 & Rad > 0)  = id_x_m1(Rad < p1 & Rad > 0);
        id_x(Rad < p2 & Rad > p1) = id_x_p1(Rad < p2 & Rad > p1);
        id_y(Rad < p3 & Rad > p2) = id_y_m1(Rad < p3 & Rad > p2);
        id_y(Rad < p4 & Rad > p3) = id_y_p1(Rad < p4 & Rad > p3);

        %% based on Chemotaxis Motion
        
        % Compute next steps
        id_x_p1 = id_x+1;
        id_x_m1 = id_x-1;
        id_y_p1 = id_y+1;
        id_y_m1 = id_y-1;
        id_x_p1(id_x_p1 > params.Nx) =  params.Nx;
        id_x_m1(id_x_m1 < 1) =  1;
        id_y_p1(id_y_p1 > params.Ny) =  params.Ny;
        id_y_m1(id_y_m1 < 1) =  1;

        here_indices = sub2ind(size(Ic), id_y, id_x);
        left_indices = sub2ind(size(Ic), id_y, id_x_m1);
        right_indices = sub2ind(size(Ic), id_y, id_x_p1);
        down_indices = sub2ind(size(Ic), id_y_m1, id_x);
        up_indices = sub2ind(size(Ic), id_y_p1, id_x);

        C_Here =  C(here_indices);
        C_Left =  C(left_indices);
        C_Right = C(right_indices);
        C_Down =  C(down_indices);
        C_Up =    C(up_indices);

        % Compute movement probabilities
        prob_left =  max(0, (1-Ic(left_indices)).*(C_Left-C_Here).*(params.C_chi/4) );
        prob_right = max(0, (1-Ic(right_indices)).*(C_Right-C_Here).*(params.C_chi/4) );
        prob_down =  max(0, (1-Ic(down_indices)).*(C_Down-C_Here).*(params.C_chi/4) );
        prob_up =    max(0, (1-Ic(up_indices)).*(C_Up-C_Here).*(params.C_chi/4) );
        prob_stay = (1-prob_left - prob_right - prob_down - prob_up);

        % Normalize probabilities (total probability of moving = 1)
        total_prob = prob_left + prob_right + prob_down + prob_up + prob_stay;
        prob_left = prob_left ./ total_prob;
        prob_right = prob_right ./ total_prob;
        prob_down = prob_down ./ total_prob;
        prob_up = prob_up ./ total_prob;
        prob_stay = prob_stay./total_prob;

        p1 = prob_left;
        p2 = prob_left + prob_right;
        p3 = prob_left + prob_right + prob_down;
        p4 = prob_left + prob_right + prob_down + prob_up;
        

        Rad=rand(params.num_walkers,1);
        % Choose a movement direction based on probabilities
        % id_x(Rad < p1 & Rad > 0)  = id_x(Rad < p1 & Rad > 0)  - 1;
        % id_x(Rad < p2 & Rad > p1) = id_x(Rad < p2 & Rad > p1) + 1;
        % id_y(Rad < p3 & Rad > p2) = id_y(Rad < p3 & Rad > p2) - 1;
        % id_y(Rad < p4 & Rad > p3) = id_y(Rad < p4 & Rad > p3) + 1;
        % 
        % id_x(id_x > params.Nx) =  params.Nx;
        % id_x(id_x < 1) =  1;
        % id_y(id_y > params.Ny) =  params.Ny;
        % id_y(id_y < 1) =  1;

        id_x(Rad < p1 & Rad > 0)  = id_x_m1(Rad < p1 & Rad > 0);
        id_x(Rad < p2 & Rad > p1) = id_x_p1(Rad < p2 & Rad > p1);
        id_y(Rad < p3 & Rad > p2) = id_y_m1(Rad < p3 & Rad > p2);
        id_y(Rad < p4 & Rad > p3) = id_y_p1(Rad < p4 & Rad > p3);

        %% Update the walker's position
        walker_positions(:,1) = params.dx*(id_x-1);
        walker_positions(:,2) = params.dy*(id_y-1);

        %% Update the activation level
        states_can_change = ( walker_activation < params.activatedAge);

        here_indices = sub2ind(size(Ic), id_y, id_x);

        i_DC_id = BoundaryDC(here_indices) ;
        near_boundary_id = find(i_DC_id > 0 & states_can_change > 0);
        near_boundary = (i_DC_id > 0);
        % Compute probabilities
        prob_up = states_can_change .* near_boundary .*       params.P_A * params.dt .* ( 1 - walker_activation./params.activatedAge);
        prob_down = states_can_change .* (1-near_boundary) .* params.P_D * params.dt .* walker_activation/params.activatedAge;
        prob_stay = 1 - prob_up - prob_down;

        % Ensure probabilities are valid
        prob_up = max(0, min(prob_up, 1));
        prob_down = max(0, min(prob_down, 1));
        prob_stay = max(0, 1 - prob_up - prob_down);

        % Draw a random event
        Rad=rand(params.num_walkers,1);
        here_indices = sub2ind(size(DCLingerTime), near_boundary_id, i_DC_id(near_boundary_id));
        DCLingerTime(here_indices) = DCLingerTime(here_indices) + params.dt;
        
        walker_activation(Rad < prob_up & Rad > 0) =  walker_activation(Rad < prob_up & Rad > 0) + 1;
        walker_activation(Rad < prob_up + prob_down & Rad > prob_up) =  walker_activation(Rad < prob_up + prob_down & Rad > prob_up) - 1;
        % walker_activation( walker_activation < 0) = 0;
        % walker_activation( walker_activation > params.activatedAge) = params.activatedAge;

    for i = 1:params.num_walkers
        id_x_i = id_x(i);
        id_y_i = id_y(i);

        U(id_y_i,id_x_i) = U(id_y_i,id_x_i) + 1;
        A(id_y_i,id_x_i) = A(id_y_i,id_x_i) + walker_activation(i)/params.activatedAge;
    end
end
