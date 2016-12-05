function observer_test(data, iters, Z, C)    
    figure
    N = size(data, 2);
    scatter3(data(1,:),data(2,:),data(3,:),50*ones(1,N),[Z; zeros(1,N)]');   
    hold on; scatter3(C(1,:),C(2,:),C(3,:),'bx');
    waitforbuttonpress;

% disp(iters);
% disp(Z);
% disp(C);
end
