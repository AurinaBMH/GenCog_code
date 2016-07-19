

WhatTypeNetwork = 'wu';
nullM = {'Randomise_weights_keeping_topology'} ;

wei_freq = 1; 
numIter = 20;
Method = 'FACT';
Parcel = {'cust250'};
for i=2:6
    Adj = squeeze(Adj_all(i,:,:));
    Adj = Adj(1:250, 1:250); 
    CheckRC(Method, Parcel, Adj, WhatTypeNetwork, nullM{1}, wei_freq, numIter);
end