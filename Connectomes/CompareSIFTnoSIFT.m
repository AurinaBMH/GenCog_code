
clear all; close all; 

cd ('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/Denoise/WithoutSIFT');
load('iFOD2_custom200ANDaseg_default_withoutdenoise.mat')
[A1] = connectomeGroupThreshold(density, 0.6, 2); 
[L1] = connectomeGroupThreshold(len, 0.6, 2); 
V1(:,1) = log(A1(:));
V1(:,2) = L1(:); 
V1 = V1(all(V1,2),:);
figure; scatter(V1(:,2), V1(:,1), '.b');title('asegNOsift');xlabel('distance'); ylabel('log(weight)');
figure; histogram(V1(:,1), 30);title('Weight distribution');

load('iFOD2_custom200ANDfirst_default_withoutdenoise.mat')
[A2] = connectomeGroupThreshold(density, 0.6, 2); 
[L2] = connectomeGroupThreshold(len, 0.6, 2); 
V2(:,1) = log(A2(:));
V2(:,2) = L2(:); 
V2 = V2(all(V2,2),:);
figure; scatter(V2(:,2), V2(:,1), '.b');title('firstNOsift');xlabel('distance'); ylabel('log(weight)');
figure; histogram(V2(:,1), 30);title('Weight distribution');

cd('/Users/Aurina/GoogleDrive/Genetics_connectome/Gen_Cog/Data/Connectomes/Denoise');
load('iFOD2_custom200ANDaseg_withoutdenoise.mat')
[A3] = connectomeGroupThreshold(density, 0.6, 2); 
[L3] = connectomeGroupThreshold(len, 0.6, 2); 
V3(:,1) = log(A3(:));
V3(:,2) = L3(:); 
V3 = V3(all(V3,2),:);
figure; scatter(V3(:,2), V3(:,1), '.b');title('asegYESsift');xlabel('distance'); ylabel('log(weight)');
figure; histogram(V3(:,1), 30);title('Weight distribution');


load('iFOD2_custom200ANDfirst_withoutdenoise.mat')
[A4] = connectomeGroupThreshold(density, 0.6, 2); 
[L4] = connectomeGroupThreshold(len, 0.6, 2); 
V4(:,1) = log(A4(:));
V4(:,2) = L4(:); 
V4 = V4(all(V4,2),:);
figure; scatter(V4(:,2), V4(:,1), '.b');title('firstYESsift');xlabel('distance'); ylabel('log(weight)');
figure; histogram(V4(:,1), 30);title('Weight distribution');

% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 1049, 895]);
% subplot(2,1,1); imagesc(log(A1)); axis square; title('asegNOsift'); subplot(2,1,2); histogram(asegNOsift, 30);
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 1049, 895]);
% subplot(2,1,1); imagesc(log(adjGrpNOdnYESsift));axis square; title('firstYESsift');subplot(2,1,2); histogram(firstYESsift, 30);
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 1049, 895]);subplot(2,1,1); imagesc(log(adjGrpYESdnNOsift));axis square; title('firstNOsift');subplot(2,1,2); histogram(firstNOsift, 30);
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 100, 1049, 895]);subplot(2,1,1); imagesc(log(adjGrpYESdnYESsift));axis square; title('asegYESsift');subplot(2,1,2); histogram(asegYESsift, 30);
