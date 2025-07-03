function [x_DC,y_DC, cluster_id] = DC_Clusters(NumDCs,R_DC,numberOfClusters,domainBoundary)
    
    cluster_id = zeros(NumDCs,1);
    cellsPerCluster = ceil(NumDCs/numberOfClusters);
    cell_ring_vec = [0 1 6 18 36 60 90 125 169 216 270];
    cell_rings = 1;
    if (cellsPerCluster > 0);        cell_rings = 2;    end
    if (cellsPerCluster > 1);        cell_rings = 3;    end
    if (cellsPerCluster > 7);        cell_rings = 4;    end
    if (cellsPerCluster > 19);       cell_rings = 5;    end
    if (cellsPerCluster > 37);       cell_rings = 6;    end
    if (cellsPerCluster > 61);       cell_rings = 7;    end
    if (cellsPerCluster > 91);       cell_rings = 8;    end
    if (cellsPerCluster > 126);      cell_rings = 9;    end
    if (cellsPerCluster > 169);      cell_rings = 10;    end
    if (cellsPerCluster > 217);      cell_rings = 11;   end
    if (cellsPerCluster > 271);      cell_rings = 12;   end
    
    N = numberOfClusters; % Number of integers
    M = NumDCs; % Desired sum
    alpha = cellsPerCluster; % Desired mean
    V = (cell_rings); % Desired variance
    ClusterSizes = random_integers_with_sum_mean_and_variance(N, M, alpha, V);
    
    % % Set the boundaries
    dist_sep = 0.25;

    x_min = domainBoundary.x_min+1;
    x_max = domainBoundary.x_max-1;
    y_min = domainBoundary.y_min+1;
    y_max = domainBoundary.y_max-1;

    x_length = x_max - x_min-2;
    y_length = y_max - y_min-2;

    boundedRegion = [0 , 0;
                     0 , 1;
                     1 , 1;
                     1 , 0];
    numb_times = 0;
    bad_top = true;
    while bad_top
        try
            numb_times = numb_times + 1;
            [x_mass, y_mass] = lloydsAlgorithm(rand(numberOfClusters,1),rand(numberOfClusters,1), boundedRegion, 50, 0);
            x_mass = x_mass*x_length + (x_min+1);
            y_mass = y_mass*y_length + (y_min+1);
            

            bad_top = false;
        catch
            if numb_times > 1000
                disp('fark')
                break
            end
        end
    end
    iter = 0;
    x_DC = zeros(NumDCs,1);
    y_DC = zeros(NumDCs,1);
    
    p1 = perms([0 -1 1]);
    p2 = [perms([2 -1 -1]); perms([2 -2 0]); perms([-2 1 1])];
    p2 = unique(p2,'rows');
    p3 = [perms([1 2 -3]); perms([-3 0 3]); perms([-2 -1 3])];
    p3 = unique(p3,'rows');
    p4 = [perms([0 4 -4]); perms([1 3 -4]); perms([2 2 -4]); perms([-2 -2 4]); perms([-3 -1 4])];
    p4 = unique(p4,'rows');
    p5 = [perms([0 5 -5]); perms([1 4 -5]); perms([2 3 -5]); perms([-1 -4 5]); perms([-2 -3 5])];
    p5 = unique(p5,'rows');
    p6 = [perms([0 6 -6]); perms([2 4 -6]); perms([1 5 -6]); perms([3 3 -6]); perms([-2 -4 6]); perms([-1 -5 6]); perms([-3 -3 6])];
    p6 = unique(p6,'rows');
    p7 = [perms([0 7 -7]); perms([3 4 -7]); perms([2 5 -7]); perms([1 6 -7]); perms([-3 -4 7]); perms([-2 -5 7]); perms([-1 -6 7]);];
    p7 = unique(p7,'rows');
    p8 = [perms([0 8 -8]); perms([1 7 -8]); perms([2 6 -8]); perms([3 5 -8]); perms([4 4 -8]); perms([-1 -7 8]); perms([-2 -6 8]); perms([-3 -5 8]); perms([-4 -4 8]);];
    p8 = unique(p8,'rows');
    p9 = [perms([0 9 -9]); perms([1 8 -9]); perms([2 7 -9]); perms([3 6 -9]); perms([4 5 -9]); perms([-1 -8 9]); perms([-2 -7 9]); perms([-3 -6 9]); perms([-4 -5 9]);];
    p9 = unique(p9,'rows');

    ps = [p1; p2; p3; p4; p5; p6; p7; p8; p9];

    sqrt3 = sqrt(3);

    
    for jj=1:numberOfClusters
        cellsPerCluster_ii = ClusterSizes(jj);
        DCs_Ass = zeros(cellsPerCluster_ii,1);
        cluster_iter = 1;
    
        if (cellsPerCluster_ii > 0)
            [x,y] = hex_to_rect(0,0,0,sqrt3);
            iter = iter+1;
            x_DC(iter) = R_DC*x; y_DC(iter) = R_DC*y;
            DCs_Ass(cluster_iter) = iter;
            cluster_id(iter) = jj;
        end
        if (cellsPerCluster_ii > 1)
            pt1 = ps(1:(cellsPerCluster_ii-1),:);
            
            if cellsPerCluster_ii > 3
                pt2 = ps(cell_ring_vec(cell_rings-1):cell_ring_vec(cell_rings),:);
                pt3 = [pt1; pt1; pt1; pt1; pt2];

                pc = choose_weighted_random_unique_numbers(pt3, (cellsPerCluster_ii-1));
            else
                pc = pt1;
            end

            for ii=1:size(pc,1)
                [x,y] = hex_to_rect(pc(ii,1),pc(ii,2),pc(ii,3),sqrt3);
                iter = iter+1;
                cluster_iter = cluster_iter+1;
                x_DC(iter) = R_DC*x; y_DC(iter) = R_DC*y;
                DCs_Ass(cluster_iter) = iter;
                cluster_id(iter) = jj;
            end
        end
        % x_DC(DCs_Ass) = x_DC(DCs_Ass) - mean(x_DC(DCs_Ass)) + x_mass(jj);
        % y_DC(DCs_Ass) = y_DC(DCs_Ass) - mean(y_DC(DCs_Ass)) + y_mass(jj);

        x_DC(DCs_Ass) = x_DC(DCs_Ass)  + x_mass(jj);
        y_DC(DCs_Ass) = y_DC(DCs_Ass)  + y_mass(jj);


    end
    
    x_DC_t = x_DC;
    y_DC_t = y_DC;
    if rand < 0.5
        y_DC = x_DC_t;
        x_DC = y_DC_t;
    end
    
end

%%

% Convert hex coordinates to rectangular
function [x,y] = hex_to_rect(u,v,w,sqrt3)
    side = 1;
    x = side*(u - v/2 - w/2);
    y = side*((v - w) * sqrt3 / 2);
end


function integers = random_integers_with_sum_mean_and_variance(N, M, alpha, V)
    % Generate N random numbers from a normal distribution
    random_numbers = randn(1, N);
    
    % Scale the random numbers to fit the desired variance
    scaled_numbers = random_numbers * sqrt(V);
    
    % Scale the mean of the scaled numbers to match alpha
    scaled_mean = mean(scaled_numbers);
    scaled_numbers = scaled_numbers + (alpha - scaled_mean);
    
    % Scale the sum of the scaled numbers to match M
    scaled_sum = sum(scaled_numbers);
    scaled_numbers = scaled_numbers * (M / scaled_sum);
    
    % Round the scaled numbers to integers
    integers = round(scaled_numbers);
    
    % Adjust the last number to make the sum exactly M
    integers(end) = integers(end) + (M - sum(integers));

    while nnz(integers > 0) < size(integers,2)
        [~,ind_max] = max(integers);
        [~,ind_min] = min(integers);
        
        integers(ind_min(1)) = integers(ind_min(1)) + 1;
        integers(ind_max(1)) = integers(ind_max(1)) - 1;
    end

end


function unique_numbers = choose_weighted_random_unique_numbers(numbers, N)
    
    % Find unique numbers and their frequencies
    [unique_numbers, ~, idx] = unique(numbers,'rows');
    frequency = accumarray(idx, 1);
    
    % Calculate weights based on frequency
    weights = frequency / sum(frequency);

    % Randomly select N unique numbers based on weights
    unique_indices = datasample(1:length(unique_numbers),N,'Replace',false,'Weights',weights);
    % unique_indices = randsample(length(unique_numbers), N, true, weights);
    unique_numbers = unique_numbers(unique_indices,:);
end



