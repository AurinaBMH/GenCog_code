
figure;
distances = [0 30 60 90 120 150]; 
for l=1:100
for i=1:length(distances)
IND1 = find(DistExpVect(:,1)>distances(i)); 
IND2 = find(DistExpVect(:,1)<=distances(i)+30); 
IND = intersect(IND1,IND2); 

ind = randi(length(IND),2000,1); 
[r(l,i,1), r(l,i,2)] = corr(DistExpVect(IND(ind),1),DistExpVect(IND(ind),2)); 
subplot(1,6,i); scatter(DistExpVect(IND(ind),1), DistExpVect(IND(ind),2),'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5);
xlabel('Distance between samples, mm'); 
ylabel('Coexpression'); 
end
end
MeanVals = squeeze(mean(r,1)); 