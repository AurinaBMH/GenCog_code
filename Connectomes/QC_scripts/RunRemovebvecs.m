
%% Run removebvecs
clear all; 

% first number - subject, second vector - volumes to be excluded
%[grad] = Removebvecs(subject, [volume1 volume2]);


[grad] = Removebvecs(38, [37]);
[grad] = Removebvecs(44, [7 24 26]);
[grad] = Removebvecs(51, [29 63]);
[grad] = Removebvecs(78, [1]);
[grad] = Removebvecs(82, [14]);
[grad] = Removebvecs(102, [37]);
[grad] = Removebvecs(116, [19 49 57]);
[grad] = Removebvecs(143, [24]);


[grad] = Removebvecs(144, [15 24 36]);
[grad] = Removebvecs(150, [64]);
[grad] = Removebvecs(166, [48]);
[grad] = Removebvecs(177, [7]);
[grad] = Removebvecs(185, [34]);
[grad] = Removebvecs(188, [32]);
[grad] = Removebvecs(191, [9]);
[grad] = Removebvecs(204, [38]);
[grad] = Removebvecs(207, [52]);
[grad] = Removebvecs(213, [62]);


[grad] = Removebvecs(243, [32]);
[grad] = Removebvecs(285, [19]);
[grad] = Removebvecs(280, [36 53]);

