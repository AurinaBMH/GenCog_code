%%%%%%%

CorrectedCoexpression = reshape(Residuals,[1285, 1285]);
ParcelCoexpression = zeros(225,225);
ROIs = CombinedExp(:,2);

W = unique(ROIs);
for i=1:length(W)
    for j=1:length(W)
        A = find(ROIs == W(i));
        B = find(ROIs == W(j));
        P = CorrectedCoexpression(A, B);
       
        ParcelCoexpression(i,j) = mean(mean(P));

    end
end

figure; imagesc(ParcelCoexpression); caxis([-1,1])
colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);
O = ParcelCoexpression(:); figure; histogram(O);