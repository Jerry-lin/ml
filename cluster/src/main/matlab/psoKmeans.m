function [idx, varargout] = psoKmeans(K, D, varargin)
% psoKmeans implementation
% function(K, D, m)
%  - K: cluster number
%  - D: data set
%  - pNum: particle number
% function(..., 'observer', fun_handle, 'maxiter', maxiter)
% observer:
%  - function handle: @(data, iters, Z, C)
%           - iter: iteration count
%           - Z: clustering labels
%           - C: clustering centers
% terminate condition:
%  - defualt maxiter: 100
% return idx - clustering of each data point
%        [idx, C] - C: clustering centers
%        [idx, C, iters] - iters: iteration count


%     default values
    default_handle = @NOP;
    defualt_threshold = 1e-13;
    defualt_maxiter = 100;
    defualt_w_max = 0.9;
    defualt_w_min = 0.1;
    default_v_max = 0.2;
    defualt_c1 = 2.0;
    defualt_c2 = 2.0;
    default_pNum = 20;
   
%     parse input parameters 
    p = inputParser;
    addRequired(p,'K', @isnumeric);
    addRequired(p,'D');
    
    addParameter(p,'v_max', default_v_max);
    addParameter(p,'pNum', default_pNum);
    addParameter(p,'observer', default_handle);
    addParameter(p,'threshold', defualt_threshold, @isnumeric);
    addParameter(p,'maxiter', defualt_maxiter, @isnumeric);
    addParameter(p,'W_max', defualt_w_max, @isnumeric);
    addParameter(p,'W_min', defualt_w_min, @isnumeric);
    addParameter(p,'c1', defualt_c1, @isnumeric);
    addParameter(p,'c2', defualt_c2, @isnumeric);
    
    parse(p, K, D, varargin{:});
   
    assert(p.Results.K > 0);
    assert(p.Results.threshold > 0);
    
    % pso initialization ----------------------------------
    [Z_min, Z_max, X, V, ...
        P_best_fitness, P_best_x, ...
        G_best_fitness, G_best_x] = initPso(p);
 
    % pso algorithm      ----------------------------------
    iter = 1;
    threshold = inf;
    while iter < p.Results.maxiter && threshold > p.Results.threshold
        
        % update w, V, X
        w = updateWeight(p.Results.W_max, p.Results.W_min, ...
            p.Results.maxiter, iter);
        V = updateVelocity(V, w, p.Results.c1, p.Results.c2, ...
            P_best_x, G_best_x, X, p.Results.v_max);
        X = updatePositon(X, V, Z_max, Z_min);
        
        % update Particle
        prev_G_best_fitness = G_best_fitness;
        updateParitcle(X, p.Results.D);
        threshold = abs(prev_G_best_fitness - G_best_fitness);
    end
 
    % k-means algorithm   --------------------------------
    [idx, C] = kmeans(p.Results.K, p.Results.D, 'maxiter', p.Results.maxiter,...
        'centers', G_best_x, 'observer', @observer_test_kmeans);
    
    varargout{1} = C;

    
%     embedded functions
    function updateParitcle(X, D)
        [nRows, gNum, particleNum] = size(X);
         
        for i = 1:particleNum
            GroupIdx = nnClustering(X(:,:,i), D);
            X(:, :, i) = generateParticle(GroupIdx, D, gNum);
            fitness = calculateFitness(X(:, :, i), GroupIdx, D);
            if fitness < P_best_fitness
                P_best_fitness = fitness;
                P_best_x(:, :, i) = X(:, :, i);
            end
        end
        
        [G_best_fitness, G_best_fitness_pos] = min(P_best_fitness);
        G_best_x = X(:, :, G_best_fitness_pos);
    end

    function GroupIdx = nnClustering(Xi, D)
        [nRows, gNum] = size(Xi);
        dNum = size(D, 2);
        S = ceil((dNum / gNum) / 2);
        GroupIdx = zeros(gNum, S);
        
        % index of temp_dis
        % 1-idx: distance to each center
        % 2-idx: current indexs of GroupIdx rows
        % 3-idx: lenghts of GroupIdx rows
        temp_dis = zeros(gNum, 3); 
        temp_dis(:, 3) = S;
        for i = 1:dNum
            for j = 1:gNum
                temp_dis(j, 1) = squareEuclideanDistance(D(:, i), Xi(:, j));
            end
            [min_dis, min_idx] = min(temp_dis(:, 1));
            
            % extend GroupIdx
            if temp_dis(min_idx, 2) == temp_dis(min_idx, 3)
                GroupIdx(:, temp_dis(min_idx, 2) + 1 : temp_dis(min_idx, 2) + S) = 0;
                temp_dis(min_idx, 3) = temp_dis(min_idx, 3) + S;
            end
            
            temp_dis(min_idx, 2) = temp_dis(min_idx, 2) + 1;
            GroupIdx(min_idx, temp_dis(min_idx, 2)) = i;
        end

    end

    function [Z_min, Z_max, X, V, ...
            P_best_fitness, P_best_x, ...
            G_best_fitness, G_best_x] = initPso(p)
        [M, N] = size(D); % M the properties of data point, N the number of data point
        assert(K < N);
        
        % calculate Z_min, Z_max
        Z_min = zeros(M, 1);
        Z_max = zeros(M, 1);
        for i = 1:M
            Z_min(i, 1) = min(p.Results.D(:, i));
            Z_max(i, 1) = max(p.Results.D(:, i));
        end
        % get the first particles
        % generate m particles, and caculating fitness
        
        Fitness = zeros(1, p.Results.pNum);
        X = zeros(M, p.Results.K, p.Results.pNum);
        V = zeros(M, p.Results.K, p.Results.pNum);
        for i = 1:p.Results.pNum
            GroupIdx = randomGroups(p.Results.K, N);
            X(:, :, i) = generateParticle(GroupIdx, p.Results.D, p.Results.K);
            V(:, :, i) = rand(1) * X(:, :, i); % 0.1 * X(:, :, i);
            V(:, :, i) = hardBoundaryVelocity(V(:, :, i), p.Results.v_max);
            Fitness(i) = calculateFitness(X(:, :, i), GroupIdx, p.Results.D);
        end
            
        P_best_fitness = Fitness;
        P_best_x = X;
  
        [G_best_fitness, G_best_fitness_pos] = min(P_best_fitness);
        G_best_x = X(:, :, G_best_fitness_pos);
    end

    function Xi = hardBoundaryPositon(Xi, Z_max, Z_min)
        [row_num, col_num] = size(Xi);
        assert(row_num == size(Z_max, 1));
        assert(row_num == size(Z_min, 1));
       
        for i = 1:col_num
            for j = 1:row_num
                if Xi(j, i) > Z_max(j)
                    Xi(j, i) = Z_max(j);
                elseif Xi(j, i) < Z_min(j)
                    Xi(j, i) = Z_min(j);
                end
            end
        end
    end

    function Vi = hardBoundaryVelocity(Vi, v_max)
        [row_num, col_num] = size(Vi);
        assert(v_max > 0);
        v_min = -v_max;
        for i = 1:col_num
            for j = 1:row_num
                if Vi(j, i) > v_max
                    Vi(j, i) = v_max;
                elseif Vi(j, i) < v_min
                    Vi(j, i) = v_min;
                end
            end
        end
    end

    function GroupIdx = randomGroups(K, dNum)
        % indexs of data that randomly classed into k groups
        % using value 0 as placeholder
        Idxs = randperm(dNum);
        i = ceil(dNum / K);
        s_mod = mod(dNum, K);
        if 0 ~= s_mod
            Idxs(dNum + 1 : K * i) = 0;
        end
        GroupIdx = reshape(Idxs, [K, i]);
    end

% get group data
    function Dat = getDatasByGroupIdx(GroupIdx, D, rowIdx)
        dRows = size(D, 1);
        gNum = size(GroupIdx, 2);
        Dat = zeros(dRows, gNum);
        t = 0;
        id = 1;
        for j = 1:gNum
            if 0 == GroupIdx(rowIdx, j)
                t = t + 1;
            else
                Dat(:, id) = D(:, GroupIdx(rowIdx, j));
                id = id + 1;
            end
        end
        if 0 ~= t
            Dat(:, gNum - t + 1 : gNum) = [];
        end
    end

% square Euclidean distance
    function dis = squareEuclideanDistance(Point1, Point2)
        dis = sum( (Point1 - Point2).^2 );
    end

    function dis = euclideanDistance(Point1, Point2)
        dis = sqrt(squareEuclideanDistance(Point1, Point2));
    end

% pso algorithm functions
    function Xi = generateParticle(GroupIdx, D, K)
        % mean of each data point in each group
        dRows = size(D, 1);
        Xi = zeros(dRows, K); % particle
        for i = 1:K
            Datas = getDatasByGroupIdx(GroupIdx, D, i);
            Xi(:, i) =  mean(Datas, 2);
        end
    end

    function fitness = calculateFitness(Xi, GroupIdx, D)
        dNum = size(D, 2);
        gNum = size(GroupIdx, 1);
        fitness = 0; % fitnesse
        for i = 1:gNum
            Datas = getDatasByGroupIdx(GroupIdx, D, i);
            count = size(Datas, 2);
            for j = 1 : count
                fitness = fitness + squareEuclideanDistance(Datas(:, j), Xi(:, i));
            end
        end
        fitness = sqrt(fitness) / dNum;
    end

    function w = updateWeight(W_max, W_min, Iter_max, iter)
        w = W_max - ( (W_max - W_min) * iter ) / Iter_max;
    end

    function V = updateVelocity(V, w, c1, c2, Pbest, Gbest, X, v_max)
        v_count = size(V, 3);
        for i = 1:v_count
            V(:, :, i) = w * V(:, :, i) + c1 * rand(1) * (Pbest(:, :, i) - X(:, :, i)) ...
                + c2 * rand(1) * (Gbest - X(:, :, i));
            V(:, :, i) = hardBoundaryVelocity(V(:, :, i), v_max);
        end
    end

    function X = updatePositon(X, V, Z_max, Z_min)
        pos_count = size(X, 3);
        for i = 1:pos_count
            X(:, :, i) = X(:, :, i) + V(:, :, i);
            X(:, :, i) = hardBoundaryPositon(X(:, :, i), Z_max, Z_min);
        end
    end
end
