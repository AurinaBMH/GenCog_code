function [Rich, Feeder, Local, MeanAll, x,prop, p1,p2,p3, t1,t2,t3] = RichDegreeBFcorrected(Rexp,Adj, ROIdist, NumBins)
% Rexp - coexpression matrix
% Adj - adjecency matrix

NumNodes = size(Adj, 1);
Degree = degrees_und(Adj);
klevel = max(Degree);
%Adjlog = logical(Adj);
%NumLinks = sum(sum(Adjlog));
ROIdist = logical(Adj).*ROIdist; 


%Y = discretize(ROIdist,linspace(min(min(nonzeros(ROIdist))), max(max(ROIdist)), NumBins+1));
% to make equally sized bins 
[A, ind] = sort(ROIdist(:), 'descend');
Info = [A,ind];
out = Info(all(Info,2),:);
out = out(200:end, :); %assign edge values accodring to row number.
sout = size(out,1); 
edgeindex = linspace(sout,1,NumBins+1); 

edges = round(out(round(edgeindex),1));

Y = discretize(ROIdist,edges);
Y(isnan(Y))=0;

Bins = unique(nonzeros(Y));

prop = zeros(klevel,length(Bins),3);

for b=1:length(Bins)
    
    Rich = zeros(klevel,1);
    Feeder = zeros(klevel,1);
    Local = zeros(klevel,1);
    MeanAll = zeros(klevel,1);
    
    Binned = zeros(size(Adj));
    Binned = (Y==b);
    %Adjlog = logical(Adj);
    Adjlog = Binned;
    %Adjlog = Adjlog.*Binned;
    %NumLinks = 1; %sum(sum(Adjlog));
    
    for k = 1:klevel
        
        isHub = Degree > k;
        mask = zeros(size(Adj));
        mask(isHub,isHub) = 1; % rich
        mask(isHub,~isHub) = 2; % feeder
        mask(~isHub,isHub) = 2;
        mask(~isHub,~isHub) = 3; % peripheral
        
        %check no zeros left
        if any(mask(:)==0)
            error('Problem making mask -- zeros remain');
        end
        % Maybe remove diagonal/lower triangle
        %%make a mask for rich, feeder and periferal connection
        
        mask = mask.*Adjlog;
        [~, p1(k),~,stats1] = ttest2(Rexp(mask==1),Rexp(mask~=1 & mask~=0), 'Vartype', 'unequal');
        t1(k) = stats1.tstat;
        [~, p2(k),~,stats2] = ttest2(Rexp(mask==2),Rexp(mask~=2 & mask~=0), 'Vartype', 'unequal');
        t2(k) = stats2.tstat;
        [~, p3(k),~,stats3] = ttest2(Rexp(mask==3),Rexp(mask~=3 & mask~=0), 'Vartype', 'unequal');
        t3(k) = stats3.tstat;
        
  
        Rich(k) = nanmean(Rexp(mask==1));
        Feeder(k) = nanmean(Rexp(mask==2));
        Local(k) = nanmean(Rexp(mask==3));
        ToMean = [Rich(k) Feeder(k) Local(k)];
        MeanAll(k) = nanmean(ToMean);
        
        R = find(mask==1); prop(k,b,1) = length(R);
        F = find(mask==2); prop(k,b,2) = length(F);
        L = find(mask==3); prop(k,b,3) = length(L);
        
        
        
        
    end
    
    
    signp1 = p1<0.05;t1sign = t1>0;
    A = [signp1; t1sign]; A=A';
    signp1 = find(A(:,1)==1 & A(:,2)==1);
    
    
    signp2 = p2<0.05;t2sign = t2>0;
    B = [signp2; t2sign]; B=B';
    signp2 = find(B(:,1)==1 & B(:,2)==1);
    
    
    signp3 = p3<0.05;t3sign = t3>0;
    C = [signp3; t3sign]; C=C';
    signp3 = find(C(:,1)==1 & C(:,2)==1);
    
    
    
    x = linspace(1,klevel,klevel);
    figure; plot(x, Rich, 'Color','r'); hold on;
    plot(x(signp1),Rich(signp1),'o','Color','r','LineWidth',3); hold on;
    plot(x, Feeder, 'Color','g'); hold on;
    plot(x(signp2),Feeder(signp2),'o','Color','g','LineWidth',3); hold on;
    plot(x, Local, 'Color','b'); hold on;
    plot(x(signp3),Local(signp3),'o','Color','b','LineWidth',3); hold on;title(sprintf('%d BIN',b));
    figure; bar(squeeze(prop(:,b,:)), 'stacked');title(sprintf('%d BIN',b));
    


    
    
    %, x, MeanAll); xlabel('Degree, k'); ylabel('Coexpression Z scored'); legend('Rich', 'Feeder', 'Local', 'MeanALL');
end
end
