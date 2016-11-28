clear
load clustering_dataset
N=size(x,2);
scatter3(x(1,:),x(2,:),x(3,:),'o')
k=waitforbuttonpress

%%%%%K-means
K=2;
%initialization
for k=1:K
    muinit(:,k)=x(:,unidrnd(N)); %pick as initial cluster center one random training sample
end
mu=muinit;
hold on
scatter3(mu(1,:),mu(2,:),mu(3,:),'rx')
k=waitforbuttonpress
Niter=10;

for nit=1:Niter
%E step
z=zeros(K,N);
for n=1:N
    for k=1:K
        a(k)=sum(abs(x(:,n)-mu(:,k)))^2;
    end
[b,c]=min(a);
z(c,n)=1;
end

%M step
for k=1:K
    mu(:,k)=(x*z(k,:)')/sum(z(k,:));
end

figure
scatter3(x(1,:),x(2,:),x(3,:),50*ones(1,N),[z; zeros(1,N)]');   
hold on; scatter3(mu(1,:),mu(2,:),mu(3,:),'rx');
k=waitforbuttonpress

end



