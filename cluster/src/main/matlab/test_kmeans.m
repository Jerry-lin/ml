clear
load 'clustering_dataset.mat'

K = 2;
data = x;
maxiter = 10;
[idx, C, iters, diff_J] = kmeans(K, data)
% [idx, C, iters, diff_J] = kmeans(K, data, 'maxiter', maxiter, 'observer', @observer_test);
% w = waitforbuttonpress;
%     if w == 0
%         disp('Button click')
%     else
%         disp('Key press')
%     end
