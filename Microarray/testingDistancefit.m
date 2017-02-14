A = DistExpVect(:,1)<120; 
DistExpVect2 = DistExpVect(A,:); 
DistExpVect = DistExpVect2; 


% remove the mean for each quantile
[xThresholds,yMeans] = BF_PlotQuantiles(DistExpVect(:,1),DistExpVect(:,2),3000,0,1); title('Coexpresion vs distance'); ylim([-0.8 1]);
Y = discretize(DistExpVect(:,1),xThresholds); 
 corrected = zeros(size(DistExpVect,1),1); 
for ii=1:size(DistExpVect,1)
    yii = Y(ii,1); 
    corrected(ii) = DistExpVect(ii,2)-yMeans(yii);
end
    
figure; scatter(DistExpVect(:,1), corrected); 



A = linspace(1,size(DistExpVect,1),size(DistExpVect,1)); 
for k=1:10
B = randsample(A,2000);
[r(k), p(k)] = corr(DistExpVect(:,1), corrected(:)); 
figure; scatter(DistExpVect(B,1), corrected(B));
end