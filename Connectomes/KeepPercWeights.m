function [AdjThr] = KeepPercWeights(Adj, perc)

A(:,1) = Adj(:);
A(:,2) = linspace(1,size(A,1), size(A,1));
% keep nonzero values with indexes in matrix
out = A(all(A,2),:);

NrWeights = round(size(out,1)*perc/100);

    [ b, ix ] = sortrows( out,1);
    B = flip(b); 

    for ii=1:NrWeights
        Weights(ii,:) = B(ii,:);
    end
    
    % create matrix full of zeros and fill in selected weights according to the index
    AdjThr = zeros(length(A),1); 
    
    for j=1:length(Weights)
        ind = Weights(j,2);
        val = Weights(j,1);
        AdjThr(ind) = val;
    end
    
    AdjThr = reshape(AdjThr, size(Adj));
    
end
