%function [lamda]=Clustering_kmwo(data, k, weight)
function [lamda]=Clustering_kmwo(data, k)

sam=data;
centerNum=k;
[samNum, samDimension] = size(sam);

% if nargin < 3
%     weight = ones(1, exampleNum);
% end

%--------------???-------------------
elementNum = zeros(1, centerNum);% ??????????
clusterElement = [];% ?????????
distance = zeros(samNum, centerNum);% ???????????????
centerC=[];mid=[];
DIM=samDimension*centerNum;
cluster_label=zeros(centerNum,samNum);%% ?????????example??????????
popsize=20;
xpop=zeros(DIM,popsize);
centerfinal=zeros(1,DIM); %% ???????k??????????
for p=1:popsize
    for i=1:centerNum
        x(i)=round(rand()*samNum);
        if (x(i)==0)
            x(i)=x(i)+1;
        end
        mid=[mid sam(x(i),:)];
    end
    xpop(:,p)=mid';%% ??????????20???
    mid=[];
end



fbest = inf;
pcen=0.1;                  % ???10%?mussels???

D1 = 1.1;                  % short-range (ssd) reference
D2 = 7.5;                  % long-range (lsd) reference
lm = 0.1;                  % scale factor for Levy walk
mu = 2.0;

% xpop = 10 * rand(DIM, popsize)-5;           % ?????popsize???????[-5,5]
fval=zeros();
maxiter = 50;

for iter = 1:maxiter       %???
    
    disp(iter);
    for ipop=1:popsize
        % fval(ipop)=feval(FUN, xpop(:,ipop)); %????????????????fval
        arxnew(ipop,:)=xpop(:,ipop)';
        for m=1:centerNum
            centerC(m,:)=arxnew(ipop,(m-1)*samDimension+1:m*samDimension);%?????????K???
        end
        centerMark(:,:,ipop)=centerC;
        
        %????
        for i = 1:(samNum)
            for j = 1:centerNum
                %%% --distance is a matrix that indicate the Euclidean distance between the example i and the center j
                distance(i, j) = sqrt(sum((sam(i, :)-centerC(j, :)).^2));
            end
        end
        
        for i = 1:samNum  %% ????? ?????????
            [value, posi] = min(distance(i, :));
            % --record the elements of each clusterElement
            elementNum(posi) = elementNum(posi) + 1;
            cluster_label(posi,elementNum(posi))=i;%% ??example????index??
            clusterElement(elementNum(posi),:,posi) = sam(i, :);
        end
        
        %??????????????????????
        % for i = 1:centerNum
        %     if(elementNum(i)==0)
        %          center(i, :)=sam(xr,:);
        %     end
        % end
        % centerMark(:,:,ipop)=center;
        
        
        % arxReset=[];
        % for i=1:centerNum
        %     arxReset=[arxReset center(i,:)];
        % end
        % xpop(:,ipop)=arxReset';
        
        
        
        %??????
        e=0;
        for i =1:centerNum
            s=0;
            for j =1:elementNum(i)
                aa=cluster_label(i,j);
                %s=s+sum((clusterElement(j,:,i)-centerC(i,:)).^2) * weight(aa);
                s=s+sum((clusterElement(j,:,i)-centerC(i,:)).^2) ;
            end
            e=e+s;
        end
        fval(ipop)=e;
        
        
        %??????
        elementNum = zeros(1, centerNum);
        clusterElement = [];
        distance = zeros(samNum, centerNum);
        
    end %% 20?????????
    
    
    [rankfval,idx] = sort(fval);%??????fval???????????????????rankfval,????fval???????idx
    
    %[fvalues] = sort(feval(FUN, xpop));     % evaluate
    for i=1:popsize
        rankPop(:,i)=xpop(:,idx(i));         %?????????rankPop
        
    end
    
    for d=1:DIM                               %????????
        tt=mod(d,samDimension);
        if tt==0
            tt=samDimension;
        end
        maxValue(d)=max(sam(:,tt));
        minValue(d)=min(sam(:,tt));
        tt=0;
        center(d)= mean(rankPop(d,1:ceil(popsize*pcen)));  % ???pcen=0.1(10%) mussels??d???(??????)????center?d?
        
        
        % Calculcate short-range (ssd) and long-range (lsd) densities for each individual
        [Apop Bpop] = ndgrid(xpop(d,:));
        %[Apop Bpop] = meshgrid(xpop(d,:));
        D = abs(Apop-Bpop);
        DD=max(D);%?????????????
        D1=D1*max(DD)/25;
        D2=D2*max(DD)/25;
        %ssd(d,:) = (sum((D < D1))-1)/ (D1*popsize);   %short-range (ssd) densities
        %lsd(d,:) = (sum((D < D2))-1)/ (D2*popsize);  %long-range (lsd) densities
        % Calculate move probability for all mussels, ????mussel?????
        %pmove(d,:) = (0.63- 1.26.*ssd(d,:) + 1.05.*lsd(d,:))>rand(1,popsize);
        %pmove(d,:) =0.63- 1.26.*ssd(d,:) + 1.05.*lsd(d,:);
        %pmove(d,:)=d*0.1*0.8*ones(1,popsize);
        % pmove(d,:)=rand(1,popsize);
        %lr(d,:)=rand(1,popsize)+1;
        lr(d,:) = lm./ (1 - rand(1,popsize)).^(1/( mu-1));        % Levy walk
        %lr(d,:) =0.1* rand.* raylrnd(1:popsize);                  % Rayleigh distribution
        % lr(d,:) = 1.* wblrnd(1,1,1,popsize);                      % Weibull distribution
    end
    
    
    % Mussel movement
    for i = 1 : popsize
        for d=1:DIM
            %if pmove(d,i)>rand()
            xpop(d,i) = xpop(d,i) +lr(d,i)*(center(d)-xpop(d,i));
            if xpop(d,i)>maxValue(d)
                xpop(d,i)=maxValue(d);
            else if xpop(d,i)<minValue(d)
                    xpop(d,i)=minValue(d);
                end
            end
            %else
            %    xpop(d,i) = xpop(d,i) +1*(rand()-0.5)*xpop(d,i);
            %xpop(d,i) = xpop(d,i);
            
            % end
        end
    end
    
    if fbest > rankfval(1)                    % keep best
        fbest = rankfval(1);
        xbest = rankPop(:,1);
    end
    %  if feval(FUN, 'fbest') < ftarget         % COCO-task achieved
    %     break;                                 % (works also for noisy functions)
    %  end
    centerfinal=rankPop(:,1)';
    
    %lr
    
end  %200??????

%--?????fljg---
fljg=zeros(1,samNum);
centers=zeros(centerNum,samDimension);
for m=1:centerNum
    centers(m,:)=centerfinal(1,(m-1)*samDimension+1:m*samDimension);% ?????????K???
end

%????
for i = 1:(samNum)
    for j = 1:centerNum
        %%% --distance is a matrix that indicate the Euclidean distance between the example i and the center j
        distance(i, j) = sqrt(sum((sam(i, :)-centers(j, :)).^2));
    end
end

for i = 1:samNum
    [value, posi] = min(distance(i, :));
    % --record the elements of each clusterElement
    fljg(i)=posi;
end

lamda=fljg;

end


