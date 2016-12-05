function [idx, varargout] = kmeans(K, data, varargin)
% k-means implementation
% function(K, data)
% function(..., 'observer', fun_handle, 'threshold', threshold, 'maxiter', maxiter)
% observer:
%  - function handle: @(data, iters, Z, C)
%           - iter: iteration count
%           - Z: clustering labels
%           - C: clustering centers
% terminate condition:
%  - defualt threshold: 1e-15
%  - defualt maxiter: 100
% return idx - clustering of each data point
%        [idx, C] - C: clustering centers
%        [idx, C, iters] - iters: iteration count
%        [idx, C, iters, diff_J] - diff_J: difference of the current and the last J

%     default values
    default_handle = @NOP;
    defualt_threshold = 1e-15;
    defualt_maxiter = 100;
   
%     parse input parameters 
    p = inputParser;
    addRequired(p,'K', @isnumeric);
    addRequired(p,'data');
    addParameter(p,'observer', default_handle);
    addParameter(p,'threshold', defualt_threshold, @isnumeric);
    addParameter(p,'maxiter', defualt_maxiter, @isnumeric);
    parse(p, K, data, varargin{:});
   
    assert(p.Results.K > 0)
   
%     size of the dataset
    [M, N] = size(p.Results.data);
    assert(p.Results.K < M * N);
    assert(p.Results.threshold > 0);
   
   
%     initialization centers
    u = initCenters();
    
%     calculate J
    prevj = inf;
    for iters = 1:p.Results.maxiter
        r = assignToCenter(u);
        u = calcCenter(r, u);
        
        j = calc_J(r, u);

%         call the observer handle
        p.Results.observer(p.Results.data, iters, r, u);
        
%         checking if the termination condition is satisfied
        diff_j = prevj - j;
        if diff_j < p.Results.threshold
            break  
        end
        prevj = j;
    end

% returns
    idx = zeros(1, N);
    for x = 1:N
        for y = 1:p.Results.K
            if r(y, x) ~= 0
                idx(x) = y;
            end
        end
    end
    varargout{1} = u;
    varargout{2} = iters;
    varargout{3} = diff_j;
    
    
%     embedded functions
   
    function c = initCenters()
       c = zeros(M, p.Results.K);
       for k = 1:p.Results.K
           c(:, k) = p.Results.data(:, unidrnd(N)); 
       end
    end

    function z = assignToCenter(u)
%         assign x to the nearest center
        z = zeros(p.Results.K, N);
        for n = 1:N
            a = zeros(1, p.Results.K);
            for k = 1:p.Results.K
%                 a(k)=sum(abs(p.Results.data(:,n) - u(:,k)))^2;
                a(k) = sum( (p.Results.data(:, n) - u(:, k)).^2 );
            end
            [b, i] = min(a);
            z(i, n) = 1;
        end
    end

    function u = calcCenter(r, u)
%         recalculate the centers
        for k = 1:p.Results.K
            u(:, k) = ( p.Results.data * r(k, :)' ) / sum(r(k,:));
        end
    end

    function j = calc_J(r, u)
%         calculate the value of distortion function J
        j = 0;
        for n = 1:N
            for k = 1:p.Results.K
                j = j + sum( r(k, n)' * (p.Results.data(:, n) - u(:, k)).^2 );
            end
        end
    end

end