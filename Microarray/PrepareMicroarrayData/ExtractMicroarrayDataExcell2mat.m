%% Author: Aurina
%% Date modified: 2016-03-22
%% Date modified: 2017-07-12
%% This script:
%   1. Loads all microarray data from excell files for each subject
%   2. Excludes custom probes;
%   3. Excludes probes with missing entrezIDs
%   4. Saves expression data, coordinates, sample structure names for all samples
%   5. Saves data for separate subjects to DataTable
%   6. Saves data for all subjects combined as variables 'MicorarrayData.mat' file
%% choose options (if you'd like to exclude braintem and cerebellum data
%points before choosing most relevant probes: ExcludeCBandBS = 1; if not ExcludeCBandBS = 0;
removeCUSTprobes = 1;
excludeCBandBS = 1;

%% load probe information (same for all subjects)
cd ('data/genes/rawData');
fprintf(1,'Loading Probes.xlsx file\n')
FileProbes = 'Probes.xlsx';
% in Probes.xlsx file probes from 1070378 down are deleted as there are no information about gene samples.
ProbeTable = readtable(FileProbes);
ProbeID = ProbeTable.probe_id;
EntrezID = ProbeTable.entrez_id;
ProbeName =  ProbeTable.probe_name;
GeneID = ProbeTable.gene_id;
GeneSymbol = ProbeTable.gene_symbol;
GeneName = ProbeTable.gene_name;
%------------------------------------------------------------------------------
%Remove probes:
%------------------------------------------------------------------------------
if removeCUSTprobes
    % Remove all CUST probes (assign NaN values for all custom probes)
    fprintf(1,'Removing CUST probes\n')
    custProbes = strfind(ProbeName, 'CUST');
    remInd = not(cellfun('isempty', custProbes));
    nCUST = sum(remInd);
    fprintf(1,'%d CUST probes removed\n', nCUST)
    ProbeName(remInd) = {NaN};
    ProbeID(remInd) = NaN;
end

% assign NaN values for all probes with missing entrezIDs
% this is the final list of probes to be used in max var calculations
ProbeID(isnan(EntrezID)) = NaN;
nEntrez = sum(isnan(EntrezID));
fprintf(1,'%d probes without entrez IDs removed\n', nEntrez)
% creat a Data cell to store the output
headerdata = {'Expression' , 'MMcoordinates', 'StructureName', 'MRIvoxCoordinates'};
headerprobe = { 'ProbeID', 'EntrezID','ProbeName', 'GeneID', 'GeneSymbol', 'GeneName'};
Data = cell(6,4);
DataProbe = cell(1,6);

%% go to each subject's directory and take the data
for subj=1:6
    fprintf(1,'Loading data for %u subject\n', subj)
    folder = sprintf('normalized_microarray_donor0%d', subj);
    cd (folder);
    %% load information specific for each subject
    FileMicroarray = 'MicroarrayExpression.csv';
    FileAnnot = 'SampleAnnot.xlsx';
    Expression = csvread(FileMicroarray);
    % rows from 57860 are removed corresponding to probes from 1070378 down as there are no information about gene samples.
    %Expression = removerows(Expression,57860:size(Expression,1));
    Expression(:,1) = [];                         % exclude probe IDs from expression matrix
    [~,~,SlabType] = xlsread(FileAnnot, 'D:D');
    [~,~, StructureName] = xlsread(FileAnnot, 'F:F');
    SlabType(1) = [];                           % remove headline
    StructureName(1) = [];                      % remove headline
    MMcoordinates = xlsread(FileAnnot, 'K:M');
    MRIvoxCoordinates = xlsread(FileAnnot, 'H:J');
    
    % keep only non custom probes with existing entrezIDs
    Expression(isnan(ProbeID),:) = [];
    
    % To exclude expression data and coordinates for braintem (BS) and cerebellum (CB)
    % exclude columns in expression and rows in coordinates if slabtype is CB or BS
    if excludeCBandBS
        fprintf('Excluding brainstem and cerebellum data\n')
        BSsamples = strfind(SlabType, 'BS');
        CBsamples = strfind(SlabType, 'CB');
        remIndBSCB = logical(not(cellfun('isempty',  BSsamples))+not(cellfun('isempty',  CBsamples)));
        nBSCB = sum(remIndBSCB);
        fprintf(1,'%d samples from brainstem and cerebellum removed\n', nBSCB)
        %remIndCB = not(cellfun('isempty',  BSsamples));
        Expression(:,remIndBSCB) = NaN;
        MMcoordinates(remIndBSCB,:) = NaN;
        MRIvoxCoordinates(remIndBSCB,:) = NaN;
        StructureName(remIndBSCB) = {NaN};
    end
    
    % keep only existing expression values
    Expression = Expression(:,all(~isnan(Expression)));
    % keep only existing coordinates
    MMcoordinates = MMcoordinates(all(~isnan(MMcoordinates),2),:); % for nan rows
    MRIvoxCoordinates = MRIvoxCoordinates(all(~isnan(MRIvoxCoordinates),2),:); % for nan rows
    % keep only existing structure names
    StructureName(cellfun(@(StructureName) any(isnan(StructureName)),StructureName)) = [];
    
    % assign output to Data cell;
    Data{subj,1} = Expression;
    Data{subj,2} = MMcoordinates;
    Data{subj,3} = StructureName;
    Data{subj,4} = MRIvoxCoordinates;
    cd ..
    %
end
%% keep only existing ProbeNames, EntrezIDs and ProbeIDs and other gene related information.
fprintf(1,'Removing irrelevant probes\n')

ProbeName(isnan(ProbeID)) = [];
EntrezID(isnan(ProbeID)) = [];
GeneID(isnan(ProbeID)) = [];
GeneSymbol(isnan(ProbeID)) = [];
GeneName(isnan(ProbeID)) = [];
ProbeID(isnan(ProbeID)) = [];

%% assign ProbeIDs, EntrezIDs and ProbeNames to Data cell.

DataProbe{1,1} = ProbeID;
DataProbe{1,2} = EntrezID;
DataProbe{1,3} = ProbeName;
DataProbe{1,4} = GeneID;
DataProbe{1,5} = GeneSymbol;
DataProbe{1,6} = GeneName;


%% make a table from all the data
DataTable = dataset({Data, headerdata{:}});
DataTableProbe = dataset({DataProbe, headerprobe{:}});

%% combine expression and coordinate values for all subjects
fprintf(1,'Combining data for all subjects\n')
Expressionall = horzcat(DataTable{1,1}, DataTable{2,1}, DataTable{3,1}, DataTable{4,1}, DataTable{5,1}, DataTable{6,1});
Coordinatesall = vertcat(DataTable{1,2}, DataTable{2,2}, DataTable{3,2}, DataTable{4,2}, DataTable{5,2}, DataTable{6,2});
StructureNamesall = vertcat(DataTable{1,3}, DataTable{2,3}, DataTable{3,3}, DataTable{4,3}, DataTable{5,3}, DataTable{6,3});
MRIvoxCoordinatesAll = vertcat(DataTable{1,4}, DataTable{2,4}, DataTable{3,4}, DataTable{4,4}, DataTable{5,4}, DataTable{6,4});

%% save relevant variables to a MicroarrayData.mat file
cd ../processedData
if ~removeCUSTprobes
    fprintf(1,'Saving data with CUST probes to the file\n')
    save('MicroarrayDataWITHCUST.mat', 'DataTable','DataTableProbe', 'Expressionall', 'Coordinatesall', 'StructureNamesall', 'MRIvoxCoordinatesAll');
    
else
    fprintf(1,'Saving data without CUST probes to the file\n')
    save('MicroarrayData.mat', 'DataTable','DataTableProbe', 'Expressionall', 'Coordinatesall', 'StructureNamesall', 'MRIvoxCoordinatesAll');
    
end