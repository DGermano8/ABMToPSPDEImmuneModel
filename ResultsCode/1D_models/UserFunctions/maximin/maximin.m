function x = maximin(n, p, varargin)

%     X = MAXIMIN(N, P, 'PARAM1', 'VAL1', ..., 'PARAM5', 'VAL5'): Generate
%     'N' uniformly-sampled points in 'P' dimensions. Valid arguments are:
% 
%     - 'iterations'        Specifies the number of iterations, per cycle.
%     - 'cycles'            Specifies the number of algorithm repetitions.
%     - 'criterion'         Specifies 'cartesian' or 'spherical' sampling.
%     - 'data'              Specifies existing data points in [M, P] form.
%     - 'initial'           Specifies a set of [N, P] points to reiterate.
% 
%     By default, these arguments are set to 10 cycles of 1000 iterations,
%     sampling points in a 'cartesian' space with no existing data points.


    % Parse inputs
    if nargin > 2
        [varargin{:}] = convertStringsToChars(varargin{:});
    end
    okargs = {'iterations', 'cycles', 'criterion', 'data', 'initial'};
    defaults = {1e3, 10, 'cartesian', [], []};
    [iter, cycles, crit, data, ini] = internal.stats.parseArgs(okargs, ...
                                      defaults, varargin{:});
    if ~isempty(ini)
        [n, p] = size(ini);
        x = ini;
        cycles = 1;
    end

    % Special cases
    if n == 1
        x = rand(1, p);
        if isempty(data)
            return
        end
    end
    if p == 1
        x = (0 : 1 / (n - 1) : 1)';
        if isempty(data)
            return
        end
    end

    % Initialise design space
    old = 0;
    best = old;
    for seed = 1 : cycles
        if isempty(ini)
            if strcmp(crit, 'cartesian')
                x = rand(n, p);
            else
                t = 2 * pi * rand(n, p - 1);
                if p == 2
                    t = [cos(t), -sin(t)];
                elseif p == 3
                    t = [sin(t(:, 1)) .* cos(t(:, 2)), sin(t(:, 1)) .* ...
                         sin(t(:, 2)), cos(t(:, 1))];
                end
                x = sqrt(rand(n, 1)) .* t;
            end
        end

        % Get all combinations of pairs
        y = [x; data];
        m = n + size(data, 1);
        t = max(1, round(m * (m - 1) / 2));
        q = zeros(t, 2);
        i = 1 : 2;
        q(1, :) = i;
        for j = 2 : t
            if i(2) < m
                q(j, 1) = i(1);
                i(2) = i(2) + 1;
                q(j, 2) = i(2);
            elseif i(1) < m - 1
                i(1 : 2) = i(1) + (1 : 2);
                q(j, 1 : 2) = i(1 : 2);
            end
        end

        % Identify adjacent points
        d = y(q(:, 1), :) - y(q(:, 2), :);
        d = sum(d .* d, 2);
        C = zeros(m);
        for j = 1 : m
            k = max(mink(d(any(q == j, 2)), 2 * p));
            k = any(q == j, 2) .* (d <= k) .* q;
            C(j, k(k > 0 & k ~= j)) = 1;
        end

        % Permutate distances
        count = 0;
        for i = 1 : iter
            t = 2 * pi * rand(n, p - 1);
            if p == 2
                t = [cos(t), sin(t)];
            else
                t = [sin(t(:, 1)) .* cos(t(:, 2)), sin(t(:, 1)) .* ...
                     sin(t(:, 2)), cos(t(:, 1))];
            end
            X = x + sqrt(rand(n, 1)) / n .* t;

            % Identify outlying points
            if strcmp(crit, 'cartesian')
                idx = ~(~any(X < 0, 2) .* ~any(X > 1, 2));
            else
                idx = any(sum(X .* X, 2) > 1, 2);
            end

            % Correct outlying points
            while any(idx, 'all')
                t = 2 * pi * rand(sum(idx), p - 1);
                if p == 2
                    t = [cos(t), sin(t)];
                else
                    t = [sin(t(:, 1)) .* cos(t(:, 2)), sin(t(:, 1)) .* ...
                         sin(t(:, 2)), cos(t(:, 1))];
                end
                X(idx, :) = x(idx, :) + sqrt(rand(sum(idx), 1)) / n .* t;
                if strcmp(crit, 'cartesian')
                    idx = ~(~any(X < 0, 2) .* ~any(X > 1, 2));
                else
                    idx = any(sum(X .* X, 2) > 1, 2);
                end
            end

            % Update best samples
            y = [x; data];
            d = zeros(n, 2);
            for j = 1 : n
                k = find(C(j, :));
                d(j, 1) = min(sum((x(j, :) - y(k, :)) .^ 2, 2));
                d(j, 2) = min(sum((X(j, :) - y(k, :)) .^ 2, 2));
            end
            idx = d(:, 2) > d(:, 1);
            x(idx, :) = X(idx, :);
            count = count + sum(idx);

            % Update connectivity
            if count >= 10 * m
                y = [x; data];
                d = y(q(:, 1), :) - y(q(:, 2), :);
                d = sum(d .* d, 2);
                for j = 1 : n
                    k = max(mink(d(any(q == j, 2)), 2 * p));
                    k = any(q == j, 2) .* (d <= k) .* q;
                    C(j, k(k > 0 & k ~= j)) = 1;
                end
                count = 0;
            end
        end

        % Update global results
        d = y(q(:, 1), :) - y(q(:, 2), :);
        d = sum(d .* d, 2);
        new = min(d);
        if new > old
            best = x;
            old = new;
        end
    end
    x = best;
end