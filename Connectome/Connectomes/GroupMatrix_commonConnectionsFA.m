
%parcelation = {'FACT_aparcaseg', 'FACT_aparcaseg', 'iFOD2_aparcaseg', 'iFOD2_aparcaseg'};

clear all; close all; 
cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/Final88/');
%Parcellation = {'cust100_FA_MD'};
Tract = {'FACT_'};
Parcellation = {'custom500ANDaseg30'};
%Parcellation = {'aparcANDaseg', 'custom200ANDaseg20', 'custom200ANDaseg', 'custom200ANDfslatlas'};

for parc = Parcellation
    
Name = strcat(Tract{1}, parc{1});
load (sprintf('%s.mat', Name));

NumSubj = size(density,2);

k=1;
density_all = zeros(NumSubj,size(Parcellation,2));
density_corrected_all = zeros(NumSubj,size(Parcellation,2));
ZeroSubjAll = cell(size(Parcellation,2),1);
low_density_ID = cell(size(Parcellation,2),1);


            
            NumNodes = size(density{1,1},1);
            
            allMatrixes = zeros(NumSubj,NumNodes, NumNodes);
%% put all matrices in a 3D matrix (subject, nodes, nodes)
        for subject = 1:NumSubj;
            %Mask = count{subject};
            % changed to count
            %density{subject}(Mask<10) = 0;
            A = density{1, subject};
            A(isnan(A)) = 0 ;
            allMatrixes(subject,:,:) = A;
        end


%% exclude low density subjects
        Adjdensity = zeros(size(allMatrixes,1),1);
        %density_corrected = zeros(NumSubj,1);

            for subj = 1: size(allMatrixes,1)
                Adjdensity(subj) = density_und(squeeze(allMatrixes(subj,:,:)));
            end
        
            low = mean(Adjdensity) - 2*std(Adjdensity);
            density_corrected = Adjdensity;
            density_corrected(Adjdensity<low) = 0;
            ZeroSubjAll{k,1} = find(density_corrected==0);
            density_corrected_all(:,k) = density_corrected;
   
            %ID=sub_id;
            allMatrixes_new = allMatrixes;
        if ~isempty(ZeroSubjAll{k,1})
                for sub=ZeroSubjAll{k}
                    allMatrixes_new(sub,:,:)=[];
                    allMatrixes=allMatrixes_new;
                    %low_density_ID{k,:} = ID{sub};
                end
        end
                density_all(:,k) = Adjdensity;
        
      %% threshold matrixes for mots common connections      
        popularity = [0 0.3 0.5 0.6 0.75 0.9];
                    if strcmp(parc{1}, 'aparcANDaseg_FA') || strcmp(parc{1}, 'aparcaseg_HCP') 
                        Adj_all = zeros(size(popularity,2),82,82);
                    elseif strcmp(parc{1}, 'custom200ANDaseg20_FA') || strcmp(parc{1}, 'cust100_HCP') 
                        Adj_all = zeros(size(popularity,2),220,220);
                    elseif strcmp(parc{1}, 'custom200ANDaseg_FA')
                        Adj_all = zeros(size(popularity,2),214,214);
                    elseif strcmp(parc{1}, 'custom200ANDfslatlas_FA')
                        Adj_all = zeros(size(popularity,2),228, 228);
                    elseif strcmp(parc{1}, 'custom500ANDaseg30') || strcmp(parc{1}, 'custom500ANDaseg30_FA')
                        Adj_all = zeros(size(popularity,2),530, 530);
                        
                    end
% saving all popularity thresholds

        q=1; % counts thresholds
        while q<size(popularity,2)

                        numMatrixes = size(allMatrixes,1);
                        numNodes = size(allMatrixes(1,:,:),2);
                        for pop = popularity;
                        %% calculating how many subjects have connections
                        Logical = zeros(numMatrixes,numNodes, numNodes);
                            for i=1:numMatrixes
                                Logical(i,:,:)=logical(allMatrixes(i,:,:)); % is the connection present (despite weight)
                                average_log=squeeze(mean(Logical, 1));      % how many % of subjects have this connection
                            end

                        Popular_connections=average_log>pop;            % calculates most popular connections
                        average_weighted=squeeze(mean(allMatrixes,1));  % calculating average of all weighted matrices       

                    Matrix_pop = zeros(numMatrixes,numNodes,numNodes);
                   % keep only popular connections in each matrix
                            for v=1:numMatrixes
                                Matrix = squeeze(allMatrixes(v,:,:));
                                Matrix_pop(v,:,:) = Matrix.*Popular_connections;
                            end
                    % calculate weights excluding zeros
                    Adj = squeeze(sum(Matrix_pop,1) ./ sum(Matrix_pop(:,:,:)~=0));
                    Adj(isnan(Adj)) = 0;
                    Adj_all(q,:,:)=Adj;
                    Adj = squeeze(Adj_all(4,:,:));
                    q=q+1;
                        end
        end
            
%cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/FinalCorrections/FA_MD/');        
save(sprintf('%s%s_CommonConnections_density.mat',Tract{1}, parc{1}),'Adj_all', 'Adj');
k=k+1;
end   
