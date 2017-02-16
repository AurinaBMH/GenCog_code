

for i=1:88

deg = degrees_und(FA{i});
Deg(i,:) = deg;
D(i,:) = zscore(deg);


end

MeanD = mean(D,1);
SDD = std(D,0,1);

MeanDeg = mean(Deg, 1);
SDDeg = std(Deg,0,1);

x = 1:1:82;
figure; bar(x,MeanDeg);
hold on; errorbar(MeanDeg,SDDeg, 'r.');
xlabel('Node','FontSize',18,'FontWeight','bold'); ylabel('Degree', 'FontSize',18,'FontWeight','bold');


figure; bar(x,MeanD);
hold on; errorbar(MeanD,SDD, 'r.');
xlabel('Node','FontSize',18,'FontWeight','bold'); ylabel('Z-scored degree', 'FontSize',18,'FontWeight','bold');




