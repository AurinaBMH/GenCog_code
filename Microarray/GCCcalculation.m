% This script is calculating gene coexpression contribution for each gene. 
% approach used in Fulcher (2015) PNAS.
%
% Aurina Arnatkeviciute, Monash University, Aug 2016.
%==========================================================================


GCC = cell(size(SelectedGenes,2),1);
for a=1:size(SelectedGenes,2)
    GCC{a} = zeros(size(SelectedGenes,1)); 
    for i=1:size(SelectedGenes,1)
        for j=i+1:size(SelectedGenes,1)
            GCC{a}(i,j) = SelectedGenes(i,a).*SelectedGenes(j,a);
            
        end
    end
    GCC{a} = GCC{a}+GCC{a}';
end


% average GCC scores according to a ROI. 

W = unique(ExpSampNormalisedAll(:,1));
ROIs = ExpSampNormalisedAll(:,1);
[sROIs, ind] = sort(ROIs);
GCCsorted = GCC(:,ind, ind);


for a=1:size(SelectedGenes,2)
    
    for i=1:length(W)
    for j=1:length(W)

        A = find(sROIs == W(i));
        B = find(sROIs == W(j));
        %for corrected
        
        P = GCCsorted(a, A, B);
        %for uncorrected
        %P = CoexpressionSorted(A, B);
        GCCroi(a,i,j) = mean(mean(P));
    end
    end
end



% on ranks
for i=1:size(SelectedGenes,2)
   GenesRanked(:,i) = floor(tiedrank(SelectedGenes(:,i)));
end


GCCrank = cell(size(SelectedGenes,2),1);
for a=1:size(SelectedGenes,2)
    GCCrank{a} = zeros(size(SelectedGenes,1));
    for i=1:size(SelectedGenes,1)
        for j=i+1:size(SelectedGenes,1)
            GCCrank{a}(i,j) = (((GenesRanked(i,a) - mean(GenesRanked(:,a)))/max(GenesRanked(:,a))))*(((GenesRanked(j,a) - mean(GenesRanked(:,a)))/max(GenesRanked(:,a))));
            
        end
    end
    GCCrank{a} = GCCrank{a}+GCCrank{a}';
end


GCCranksorted = GCC(:,ind, ind);
for a=1:size(SelectedGenes,2)
    
    for i=1:length(W)
    for j=1:length(W)

        A = find(sROIs == W(i));
        B = find(sROIs == W(j));
        %for corrected
        
        P = GCCranksorted(a, A, B);
        %for uncorrected
        %P = CoexpressionSorted(A, B);
        GCCrankroi(a,i,j) = mean(mean(P));
    end
    end
end
