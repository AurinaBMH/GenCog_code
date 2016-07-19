%% script to run Ruch club for a separate case
% cortex (is chosen FACTcust100_60_2hem_wu_cortex.mat file)
% rich club based on strength if WhatTypeNetwork = 'wus';
% addpath(genpath('/gpfs/M2Scratch/Monash076/aurina/BCT/'))
% 
 Method = 'iFOD2';
 Parcel = 'cust100';

%load(sprintf('%s%s_60_2hem_wu.mat', Method, Parcel));

deg = degrees_und(Adj);
strength = sum(Adj);
smax = floor(max(strength));
smin = min(strength);
numIter = 0; % 0 will keep topology the same
numRepeats = 150;
WhatTypeNetwork = 'wus';
wei_freq = 1; % closer to 1 will keep strength sequence more similar
%Hemisphere = '2hemb;
whatNullModel = 'null_model_und_sign';
% whatNullModel = 'shuffleWeights'; - preserves topology and node strength (not weights); 
% whatNullModel = 'Randomise_weights_keeping_topology; preserves topology and weight distribution, but not the strength of each node.
% whatNullModel = 'randmio_und'; preserves degree and weight distribution; changes topology. (to be used on binary networks)

        if strcmp(WhatTypeNetwork, 'bu') || strcmp(WhatTypeNetwork, 'bu2h')
            Adj = logical(Adj);
        end

[PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm_AS(Adj,smax,numIter,numRepeats,WhatTypeNetwork, whatNullModel, wei_freq);

               svector = linspace(smin, smax, 50);
                % Compute p-values
                pValues = zeros(smax,1);
                for i = 1:50
                        pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
                end
                % Significant phi
                isSig = (pValues <= 0.05);
                PhiNormMean = zeros(size(PhiTrue));
                for i = 1:size(PhiTrue,2)
                    PhiNormMean(i) = PhiTrue(i)/mean(PhiRand(:,i));
                end
                % plot the graphs
                figure; 
                subplot(2,1,1); plot(PhiNormMean, '-','Color','r','LineWidth',2);  %xlim([min(strength) max(strength)+2]); % ylim([0.8 max(PhiNormMean)+0.05]); ...
                legend('Phi normalised', 'Location','NorthWest'); ...
                hold on;
                plot(find(isSig),PhiNormMean(isSig),'o','Color','r','LineWidth',3); %xlim([min(strength) max(strength)+2]); % ylim([0.8 max(PhiNormMean)+0.05]);
                    title (sprintf('Normalised rich club \n %s %s', Method, Parcel)); xlabel('Strength'); ylabel('Phi (norm)'); set(gca,'FontSize',12,'fontWeight','bold');
                subplot(2,1,2); hist(strength, 20); title ('Strength distribution'); xlim([min(strength) max(strength)+2]); xlabel('Strength'); set(gca,'FontSize',12,'fontWeight','bold');
                
                
                
             