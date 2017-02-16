clear all; close all; 
%% load one file to have density for stating variables
% check if you want to cut accoring to rank or weight value!!
cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/FinalCorrections/');
 Parcellation = {'aparcaseg'};
 Tract = {'FACT_'};
% 
Name = strcat(Tract{1}, Parcellation{1});
load (sprintf('%s.mat', Name));

% choose parcellation
% fot cust100 and cust250 do not remove count<10 links

NumSubjects = size(density,2); 
NumParcel = size(Parcellation,2);

Density_adj = zeros(NumSubjects,1);
low_density_ID = cell(NumSubjects,1);
ZeroSubjAll = cell(NumSubjects,1);
density_corrected_all = zeros(NumSubjects,NumParcel);
density_all = zeros(NumSubjects,NumParcel);

k=1;
for parc = Parcellation

if strcmp(parc{1}, 'aparcaseg')
                        NumNodes = 82; 

                        
elseif strcmp(parc{1}, 'cust100')
                        NumNodes = 220; 
                        
elseif strcmp(parc{1}, 'cust250')
                        NumNodes = 530; 
                        
                                      
end


allMatrices = zeros(NumSubjects,NumNodes, NumNodes);
allMatricesCount = zeros(NumSubjects,NumNodes, NumNodes);

        for subject = 1:NumSubjects;
            Mask = count{subject};
            density{subject}(Mask<10) = 0;
            allMatrices(subject,:,:) = density{1,subject};
            allMatricesCount(subject,:,:) = count{1, subject};
        end


%% exclude low density subjects

        Adjdensity = zeros(size(allMatrices,1),1);
        % density_corrected = zeros(NumSubj,1);
            for subj = 1: size(allMatrices,1)
                Adjdensity(subj) = density_und(squeeze(allMatrices(subj,:,:)));
            end
        
            low = mean(Adjdensity) - 2*std(Adjdensity);
            density_corrected = Adjdensity;
            density_corrected(Adjdensity<low) = 0;
            ZeroSubjAll{k,1} = find(density_corrected==0);
            density_corrected_all(:,k) = density_corrected;
   
            ID=sub_id;
            allMatrices_new = allMatrices;
            if ~isempty(ZeroSubjAll{k,1})
                    for sub=ZeroSubjAll{k}
                        allMatrices_new(sub,:,:)=[];
                        allMatrices=allMatrices_new;
                        low_density_ID{k,:} = ID{sub};
                    end
            end
     density_all(:,k) = Adjdensity;
     Nonzero = sum(density_corrected>0); 
     
     
     Adj_exact = zeros(Nonzero, NumNodes,NumNodes);
     Nthresholds = zeros(Nonzero,1);
     Wthresholds = zeros(Nonzero,1);
     
     % for subjects with normal density, exclude links with weight less than
     % 10 streamlines in count. 
     for subj = 1:Nonzero
        Adj = squeeze(allMatrices(subj,:,:));
        %% exclude links with count <10
%         AdjMask = squeeze(allMatricesCount(subj,:,:));
%         AdjMask(AdjMask<10)=0;
%         AdjMask2 = logical(AdjMask);
%         Adj = Adj.*AdjMask2; 
%         %assign a matrix with removed count<10 values back to a list of all
%         %matrices.
%         allMatrices(subj,:,:) = Adj; 
        
        Weights = triu(Adj,1);
%         Counts = triu(AdjMask,1);
        WeightsSorted = sort(nonzeros(Weights));

         %% find the threshold when components start to separate
         n = 1;
         adj_temp = Adj; 
         [comps,comp_sizes] = get_components(adj_temp);

             while size(comp_sizes,2) == 1

                        indA = find(Adj > WeightsSorted(n));
                        adj_temp = zeros(size(Adj));
                        adj_temp(indA) = Adj(indA);
                        [comps,comp_sizes] = get_components(adj_temp);
                        n=n+1;

             end

        %% get the exact matrix   (-2: additional +1 was added after finding the right threshold); 
            if n==2 || n==1

                Adj_exact(subj,:,:) = Adj;
                Nthresholds(subj,1) = 0; % exclude a matrix if no links can be removed
                
            else
           
                    indexact = find(Adj > WeightsSorted(n-2));
                    adj_exact = zeros(size(Adj));
                    adj_exact(indexact) = Adj(indexact);
                    [comps_exact,comp_sizes_exact] = get_components(adj_exact);

                    Adj_exact(subj,:,:) = adj_exact;
                    Nthresholds(subj,1) = n-2;
                    Wthresholds(subj,1) = WeightsSorted(n-2);
            end
subj
     end
        
%% finding a threshold for all subjects (lowest threshold of all existing, so all matrices will be still fully connected)
    
    index = find(Nthresholds==0);
    Adj_exactFin = zeros(Nonzero,NumNodes, NumNodes);

    MinW = min(Wthresholds(Wthresholds>0));
    %MinN = min(Nthresholds(Nthresholds>0));
    %[MinInd] = find(Nthresholds == MinN);
    
    for sub = 1:Nonzero
        BadN = ismember(sub,index);
        if BadN == 1
            Adj_exactFin(sub,:,:) = 0; % exclude matrix, if threshold is zero. 
            %(as it is all 0, it will not be included in the averages ar averaging is performed only on existing connecions
        else

         Adj = squeeze(allMatrices(sub,:,:));
         Weights = triu(Adj,1);
         WeightsSorted = sort(nonzeros(Weights));


        indexact2 = find(Adj >= MinW );
        adj_exact2 = zeros(size(Adj));
        adj_exact2(indexact2) = Adj(indexact2);
        [comps_exact,comp_sizes_exact] = get_components(adj_exact2);
                
        Adj_exactFin(sub,:,:) = adj_exact2;
        Density_adj(sub,1) = density_und(adj_exact2);
        %figure; subplot(1,2,1); imagesc((log(squeeze(Adj_exact(sub,:,:)))));  subplot(1,2,2); imagesc(log((adj_exact2))) 
            
        end
    end
    numMatrices = size(Adj_exactFin,1);
    numNodes = size(Adj_exactFin(1,:,:),2);
    List = zeros(numMatrices,1);
    n=1;
    for i=1:numMatrices
        A = squeeze(Adj_exactFin(i,:,:));
        WeightSum = sum(sum(A));
        if WeightSum == 0
            List(n) = i;
        end
        n=n+1;
    end
ListToRemove = nonzeros(List);
Adj_exactFin(ListToRemove,:,:) = []; 

    %% make average matrix
    
    GroupAdj = squeeze(sum(Adj_exactFin,1) ./ sum(Adj_exactFin(:,:,:)~=0));
    GroupAdj(isnan(GroupAdj)) = 0;
    k=k+1;
     save([sprintf('%s%s_ConnectedComponents_group', Tract{1}, parc{1})], 'GroupAdj', 'allMatrices', 'Adj_exactFin', 'Density_adj', 'density_corrected', 'Nthresholds', 'Wthresholds');
end

