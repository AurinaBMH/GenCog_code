

math = density; 
[Adj] = connectomeGroupThreshold(math, 0.6); 
Adj(isnan(Adj))=0;
% for t=1:88
%     D{1,t} = pdist2(Coordinates{t}, Coordinates{t});
% end
[dist] = connectomeGroupThreshold(len, 0.6); 
dist = dist.*logical(Adj);




d = degrees_und(Adj);
% dist = pdist2(Coordinates, Coordinates);
% dist = dist.*logical(Adj);

% connection length
for j=1:size(dist,1)
    k=dist(j,:);
    D1(j) = mean(k(k~=0));
end

du = unique(d);

for i=1:length(du)
    
    A = find(d<=du(i));
    Meandist(i) = mean(D1(A));
end
figure; scatter(du,Meandist); hold on; plot(du, Meandist);xlabel('degree<=k');ylabel('Average connection distance');


% weight

for j=1:size(Adj,1)
    k=Adj(j,:);
    W1(j) = mean(k(k~=0));
end

for i=1:length(du)
    
    C = find(d<=du(i));
    Meanweight(i) = mean(W1(C));
end
figure; scatter(du,Meanweight); hold on; plot(du, Meanweight);xlabel('degree<=k');ylabel('Average connection weight');
