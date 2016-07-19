%checking RC
clear all; close all; 

cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/Final88/FA/');

Tract = {'FACT_'};
Parcellation = {'aparcANDaseg_FA', 'custom200ANDaseg20_FA', 'custom200ANDaseg_FA', 'custom200ANDfslatlas_FA'};

WhatTypeNetwork = 'bu';

if strcmp(WhatTypeNetwork, 'wu')
    whatNullModel = {'randmio_und','Randomise_weights_keeping_topology'} ;
elseif strcmp(WhatTypeNetwork, 'bu')
    whatNullModel = {'randmio_und'} ;
end
% whatNullModel = 'shuffleWeights'; - preserves topology and node strength (not weights); 
% whatNullModel = 'Randomise_weights_keeping_topology; preserves topology and weight distribution, but not the strength of each node.
% whatNullModel = 'randmio_und'; preserves degree and weight distribution; changes topology. (to be used on binary networks)
% whatNullModel = 'null_model_und_sign' ;    

for parc = Parcellation
    Name = strcat(Tract{1}, parc{1});
    load (sprintf('%s_CommonConnections_FA.mat', Name));
    for nullM = whatNullModel
        % choose options for RC
        Method = 'FACT';
        Parcel = parc;

        wei_freq = 1; 
        numIter = 20; %(0 if do not want to change topology, the higher the number, more different topology); 
       
        CheckRC(Method, Parcel, Adj, WhatTypeNetwork, nullM{1}, wei_freq, numIter);

        % get cortical matrix only
        AdjCort = GetCortexMatrix(Adj, Parcel);
        % calculate RC again
        CheckRC(Method, Parcel, AdjCort, WhatTypeNetwork, nullM{1}, wei_freq, numIter);
    end
end

