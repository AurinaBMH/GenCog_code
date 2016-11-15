% Nr = {'1004_8237'; '1006_8265'; '1007_8285'; '1010_8314'; '1012_8376'; '1014_8388'; '1018_8419'; '1020_8435'; ...
%     '1022_8457';'1024_8465'; '1026_8475'; '1029_8491'; '1032_8519'; '1034_8528'; '1036_8531'; '1038_8543'; ...
%     '1042_8596'; '1044_9033';'1046_9059'; '1048_9339'; '2008_8287'; '2009_8307'; '2013_8379'; '2015_8394'; ...
%     '2017_8413';'2019_8421';'2021_8455';'2023_8464';'2025_8471';'2027_8481';'2028_8489'; '2030_8499';'2031_8508'; ...
%     '2035_8530';'2037_8600';'2039_8544'; '2041_8589';'2043_9034';'2045_9049';'2047_9155'};

cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/MCDESPOT-MATRICES');
load('t1subjID.mat');

cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/MCDESPOT-MATRICES/T2/time1');

math = cell(1,length(Nr));
for i=1:length(Nr)
    load(sprintf('%s_HARDI_MD_C_trafo_tracts_dRL_FOD_interp_4D_new_6_END.mat',Nr{i}));
    math{1,i} = CM;
end
cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/MCDESPOT-MATRICES');
save('T2t1.mat','math'); 