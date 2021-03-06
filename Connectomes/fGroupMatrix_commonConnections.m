function Adj_all = fGroupMatrix_commonConnections(Weights, Tract,parc)

for l = 1:2
        weight = Weights{l};
        NumSubj = size(weight,2);

        % define variables
        ZeroSubjAll = cell(1,1);
        %low_density_ID = cell(size(Parcellation,2),1); 
        NumNodes = size(weight{1,1},1);
        allMatrixes = zeros(NumSubj,NumNodes, NumNodes);
        
        %% put all matrices in a 3D matrix (subject, nodes, nodes)
        for subject = 1:NumSubj;
            %Mask = count{subject};
            % remove connections with weight <10 (works for anatomical
            % parcellation only)
            %density{subject}(Mask<10) = 0;
            A = weight{1, subject};
            %replace NaN values with 0 (in the last batch there were come
            %NaN connections)
            A(isnan(A)) = 0 ;
            allMatrixes(subject,:,:) = A;
        end
        
        
        %% exclude low density subjects (density<mean - 2SD)
        Adjdensity = zeros(size(allMatrixes,1),1);
       
        for subj = 1: size(allMatrixes,1)
            % calculate density for each subject
            Adjdensity(subj) = density_und(squeeze(allMatrixes(subj,:,:)));
        end
        % define a threshold
        low = mean(Adjdensity) - 2*std(Adjdensity);
        density_corrected = Adjdensity;
        density_corrected(Adjdensity<low) = 0;
        ZeroSubjAll{1,1} = find(density_corrected==0);

        
        %remove low density matrices;
        allMatrixes_new = allMatrixes;
        if ~isempty(ZeroSubjAll{1,1})
            for sub=ZeroSubjAll{1}
                allMatrixes_new(sub,:,:)=[];
                allMatrixes=allMatrixes_new;
                %low_density_ID{k,:} = ID{sub};
            end
        end

        
        %% threshold matrixes for most common connections
        popularity = [0 0.3 0.5 0.6 0.75 0.9];
        
        Adj_all = zeros(size(popularity,2), NumNodes, NumNodes);
        % saving all thresholds
        
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
                q=q+1;
            end
        end
        
        if l==1
            weightName = {'count'};
        elseif l==2
            weightName = {'density'};
        end
        save(sprintf('%s%s_CommonConnections_%s.mat',Tract{1}, parc{1}, weightName{1}), 'Adj_all');

end
end
