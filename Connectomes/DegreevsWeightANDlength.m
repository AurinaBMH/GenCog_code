
% choose what weight measure to use;
math =FA;
GrpThr = 0; 

[Adj, ~, ~, prop] = connectomeGroupThreshold(math, GrpThr);
%Adj = Adj(1:100,1:100); 
Adj(isnan(Adj))=0;
d = degrees_und(Adj);

if size(math{1,1},1)==30 % Karen's data
    
    dist = pdist2(Coordinates, Coordinates);
    dist = dist.*logical(Adj);
    
elseif size(math,2)==88 % FA
    
    for t=1:size(math,2)
        D{1,t} = pdist2(coords{t}, coords{t});
    end
    dist = connectomeGroupThreshold(D, GrpThr);
    %dist = dist.*logical(Adj);
    %dist = dist(1:100,1:100).*logical(Adj);
    
else %count and density
    % input - distance for existing links
    [dist] = connectomeGroupThreshold(len, GrpThr);
    
end


% connection length
for j=1:size(dist,1)
    
    k=dist(j,:);
    D1(j) = mean(k(k~=0));
    
end

du = unique(d);

for i=1:length(du)
    
    A = find(d>du(i));
    Oth = find(d<=du(i));
    new = D1(A);
    other = D1(Oth);
    [h,p(i)]=ttest2(new, other, 'Vartype', 'unequal');
    Meandist(i) = mean(D1(A));
end
signp = p<0.05;
Sign = find(signp);

FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1049, 895]);


subplot(3,1,1); histogram(d,30); title ('Degree distribution'); ylabel('Frequency'); xlabel('degree'); axis([0 max(d) 0 inf]);
subplot(3,1,2);plot(du, Meandist, 'Color','r','LineWidth',2 );xlabel('degree>k');ylabel('Average connection distance'); hold on; 
plot(du(Sign),Meandist(Sign),'o','Color','r','LineWidth',3);axis([0 max(d) min(Meandist)*0.99 max(Meandist)*1.01]);
    

% weight

for j=1:size(Adj,1)
    k=Adj(j,:);
    W1(j) = mean(k(k~=0));
end

for i=1:length(du)
    
    C = find(d>du(i));
    Othw = find(d<=du(i));
    
    neww = D1(C);
    otherw = D1(Othw);
    [h,pw(i)]=ttest2(neww, otherw, 'Vartype', 'unequal');
    
    Meanweight(i) = mean(W1(C));
end

signp = pw<0.05;
Sign = find(signp);
subplot(3,1,3);plot(du, Meanweight, 'Color','r','LineWidth',2 );xlabel('degree>k');ylabel('Average connection weight'); hold on; 
plot(du(Sign),Meanweight(Sign),'o','Color','r','LineWidth',3); axis([0 max(d) min(Meanweight)*0.99 max(Meanweight)*1.01]);


Matrix = fRandom2AnatomicalParcellations_annot('custom200', 'ANDaseg', Adj); 
figure; imagesc((Matrix)); colorbar;
vect = weightDistance (dist, math, GrpThr, 0);



Vect(:,1) = dist(:); Vect(:,2) = prop(:);
%Vect = Vect(all(Vect,2),:);  %rows
figure; scatter(Vect(:,1), Vect(:,2)); xlabel('Distance'); ylabel('Connection probability');
BF_PlotQuantiles(Vect(:,1), Vect(:,2),30,0,1); xlabel('Distance'); ylabel('Connection probability');