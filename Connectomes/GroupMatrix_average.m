

cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Connectomes/FinalCorrections/');
parc = {'aparcaseg'};
Tract = {'iFOD2_'};

Name = strcat(Tract{1}, parc{1});
load (sprintf('%s.mat', Name));

if strcmp(parc{1}, 'aparcaseg')
                        NumNodes = 82; 
  
elseif strcmp(parc{1}, 'cust100')
                        NumNodes = 220; 
                        
elseif strcmp(parc{1}, 'cust250')
                        NumNodes = 530;    
                                      
end

NumSubj = size(FA,2);
allMatrixes = zeros(NumSubj, NumNodes, NumNodes);
    
    for subject = 1:NumSubj;
            %Mask = count{subject};
            %changed to count
            %count{subject}(Mask<10) = 0;
            allMatrixes(subject,:,:) = FA{1,subject};
    end
    
    Adj = squeeze(mean(allMatrixes,1));
    
    
    %save([sprintf('%s%s_Averagegroup', Tract{1}, parc{1})], 'Adj');