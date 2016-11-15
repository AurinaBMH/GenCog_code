popularity = 0.6;
allMatrices = Adj_exactFin; 
numMatrices = size(allMatrices,1);
numNodes = size(allMatrices(1,:,:),2);

                        for pop = popularity;
                        %% calculating how many subjects have connections
                        Logical = zeros(numMatrices,numNodes, numNodes);
                            for i=1:numMatrices
                                Logical(i,:,:)=logical(allMatrices(i,:,:)); % is the connection present (despite weight)
                                average_log=squeeze(mean(Logical, 1));      % how many % of subjects have this connection
                            end

                        Popular_connections=average_log>=pop;            % calculates most popular connections

                    Matrix_pop = zeros(numMatrices,numNodes,numNodes);
                   % keep only popular connections in each matrix
                            for v=1:numMatrices
                                Matrix = squeeze(allMatrices(v,:,:));
                                Matrix_pop(v,:,:) = Matrix.*Popular_connections;
                            end
                    % calculate weights excluding zeros
                    Adjgr = squeeze(sum(Matrix_pop,1) ./ sum(Matrix_pop(:,:,:)~=0));
                    Adjgr(isnan(Adjgr)) = 0;
                        end