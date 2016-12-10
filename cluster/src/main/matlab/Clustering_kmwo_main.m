clear
% data = [10.1,10.1;
%     -10.2,-10.2;
%     10.3,10.3;
%     9.9,9.9;
%     9.8,9.8;
%     -11,-11;
%     -9,-9;
%     -8,-8;
%     12,12;
%     10,10;];
% k =2;

clear
load 'clustering_dataset.mat'

data = x';
k = 2;
Clustering_kmwo(data,k)


