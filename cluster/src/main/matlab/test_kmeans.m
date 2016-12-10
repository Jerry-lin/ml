clear
load 'clustering_dataset.mat'

K = 2;
data = x;
maxiter = 10;
% [idx, C, iters, diff_J] = kmeans(K, data)
% [idx, C] = kmeans(K, data)
% [idx, C, iters, diff_J] = kmeans(K, data, 'maxiter', maxiter, 'observer', @observer_test_kmeans);
% w = waitforbuttonpress;
%     if w == 0
%         disp('Button click')
%     else
%         disp('Key press')
%     end

[idx, C] = psoKmeans(K, data)
