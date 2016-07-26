function [AdjThr] = KeepPercWeights(Adj, perc)

% function [AdjThr] = KeepPercWeights(Adj, perc)

% This function will apply a threshold to a matrix keeping a set percent of strongest nonzero weights.
%
% -------
% INPUTS:
% -------
%
% Adj  - A matrix that needs to be thresholded
%
% perc - Threshold (from small nonzero value to 100%) stating where to threshold the matrix
%
%
% -------
% OUTPUTS:
% -------
%
% AdjThr   - a thresholded matrix
% Aurina Arnatkeviciute, Monash University, July 2016.


% ----------------------------------------------------------------------------------------------------------------
% create an array of all reights in the matrix and give each value an original index value
A(:,1) = Adj(:);
A(:,2) = linspace(1,size(A,1), size(A,1));

% keep nonzero values with indexes in matrix
out = A(all(A,2),:);

% set a number of weights to retain
NrWeights = round(size(out,1)*perc/100);

[ b, ix ] = sortrows( out,1);
% descending order
B = flip(b);
% get a list of weights with indexes
for ii=1:NrWeights
    Weights(ii,:) = B(ii,:);
end

% create matrix full of zeros and fill in selected weights according to the index
AdjThr = zeros(length(A),1);

% assign weights according to indexes
for j=1:length(Weights)
    ind = Weights(j,2);
    val = Weights(j,1);
    AdjThr(ind) = val;
end

% reshape an array to a matrix
AdjThr = reshape(AdjThr, size(Adj));

end
