for i=1:15 L(i,:,:) = len{i}; end

Le = squeeze(sum(L,1) ./ sum(L(:,:,:)~=0));
A = logical(Adj);
Le = squeeze(Le);
B = A.*Le;
figure; imagesc(B);
Bl = B(:);
figure; histogram(nonzeros(Bl), 20);