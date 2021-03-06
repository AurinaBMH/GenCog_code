
%% script to run Ruch club for a separate case
% cortex (is chosen FACTcust100_60_2hem_wu_cortex.mat file)
% rich club based on strength if WhatTypeNetwork = 'wus';
% addpath(genpath('/gpfs/M2Scratch/Monash076/aurina/BCT/'))
function [] = CheckRC(Method, Parcel, Adj, WhatTypeNetwork, whatNullModel, wei_freq, numIter)


deg = degrees_und(Adj);
kmax = max(deg);

numRepeats = 100;

%WhatTypeNetwork = 'bu';
%Hemisphere = '2hemb;
%whatNullModel = 'randmio_und';
% whatNullModel = 'shuffleWeights'; - preserves topology and node strength (not weights); 
% whatNullModel = 'Randomise_weights_keeping_topology; preserves topology and weight distribution, but not the strength of each node.
% whatNullModel = 'randmio_und'; preserves degree and weight distribution; changes topology. (to be used on binary networks)

        if strcmp(WhatTypeNetwork, 'bu') || strcmp(WhatTypeNetwork, 'bu2h')
            Adj = logical(Adj);
        end

[PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm_A(Adj,kmax,numIter,numRepeats,WhatTypeNetwork, whatNullModel, wei_freq);
                
                % Compute p-values
                pValues = zeros(kmax,1);
                for i = 1:kmax
                        pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
                end
                % Significant phi
                isSig = (pValues <= 0.05);
                PhiNormMean = zeros(size(PhiTrue));
                for i = 1:size(PhiTrue,2)
                    PhiNormMean(i) = PhiTrue(i)/mean(PhiRand(:,i));
                end
                % plot the graphs
                if size(Adj,1) == 200 || size(Adj,1) == 68 || size(Adj,1) == 500
                    Part = 'Cortex';
                else 
                    Part = 'Whole Brain';
                end
                
                if strcmp(whatNullModel, 'randmio_und') && strcmp(WhatTypeNetwork, 'wu')
                    whatNullModel2 = 'WeightsANDtopology';
                elseif strcmp(whatNullModel, 'randmio_und') && strcmp(WhatTypeNetwork, 'bu')
                    whatNullModel2 = 'Topology';
                elseif strcmp(whatNullModel, 'Randomise_weights_keeping_topology')
                    whatNullModel2 = 'Weights';
                end
                    
                
                figure; 
                subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  xlim([min(deg) max(deg)+2]); % ylim([0.8 max(PhiNormMean)+0.05]); ...
                legend('Phi normalised', 'Location','NorthWest'); ...
                hold on;
                plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); xlim([min(deg) max(deg)+2]); % ylim([0.8 max(PhiNormMean)+0.05]);
                    title (sprintf('Normalised rich club \n %s %s %s %s %s', Part, Method, Parcel{1}, WhatTypeNetwork, whatNullModel2)); xlabel('Degree'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
                subplot(2,1,2); hist(deg, 20); title ('Degree distribution'); xlim([min(deg) max(deg)+2]); xlabel('Degree'); set(gca,'FontSize',12,'fontWeight','bold');
end 
                
                
             