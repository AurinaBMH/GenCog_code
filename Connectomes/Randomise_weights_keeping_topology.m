% load original weighted matrix

function [adj_new, numRewirings] = Randomise_weights_keeping_topology (Adj, numIter)

% convert original weighted matrix to logical to keep the place of existing
% connections. 
adj_orig_bin = logical(Adj);
% take all original weights and randomise them
Weights = Adj(:);
Weights(Weights==0) = [];
Weights_rand = Weights(randperm(length(Weights)));

numRewirings = 0;
%define, where the connections should be according to an old topology
adj_new = zeros(size(adj_orig_bin));

% assign randomised weights to the network.
l=1;
    for i = 1:length(adj_orig_bin)
        for j = 1:length(adj_orig_bin)
                if adj_orig_bin(i,j)==1 
                adj_new(i,j) = Weights_rand(l);
                 l=l+1;  
                end 

        end
    end
end