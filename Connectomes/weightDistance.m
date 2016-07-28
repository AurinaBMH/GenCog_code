function vect = weightDistance (coord, math, grpThr, denSD, Names)

% function vect = weightDistance (coord, math, grpThr)

% This function will calculate and plot weight-distance relationship in a
% group connectome
%
% -------
% INPUTS:
% -------
%
% coord  - a list of ROI coordinates or a distance matrix
%
% mats  - an array containing K cells each with an N*N matrix, where K =
%         number of subjects and N = number of nodes. Could optionally be a
%         3d matrix of dimensions N*N*K. 
%
% grpThr - Group-level threshold. This is the proportion of subjects who
%          must have a non-zero value for a given edge to be retained. Must
%          be between zero and one.
% denSD  - Some subjects may have a low connection density and bias the
%          results. This input allows one to exclude those subjects when
%          considering the proportion threshold, grpThr. This value should
%          represent the number of SDs above/below the mean for a person to
%          be excluded from consideration. e.g., denSD = 2 means that any
%          subject with a connectiont density +/- 2SD from mean will be
%          excluded. Note: the subject is only excluded from computing the
%          proportion of subjects with a non-zero edge value. The data from
%          the subject will still be retained and the threshold applied to
%          this subject in the outputs. To use the whole sample, set denSD
%          = 0. This is the default.
% Names - ROI labels in matlab format. 
%
% -------
% OUTPUTS:
% -------
% vect   - variable with 2 columns (distance and weight)
% scatterplot
% Aurina Arnatkeviciute, Monash University, July 2016.
%==========================================================================

%calculate pairwise distance if coordinate list is imported
if size(coord,2) == 3
    dist = pdist2(coord, coord);
else
    dist = coord;
end


% get group matrix
[adjGrp] = connectomeGroupThreshold(math, grpThr, denSD) ;
NumNodes = size(math{1,1},1);
vect(:,1) = dist(:);
vect(:,2) = adjGrp(:);


% remove zero weights
vect = vect(all(vect,2),:);


% plot matrix and distance, weight relationship
figure;
imagesc(adjGrp); colorbar;
set(gca);
ax = gca; 
ax.YLim = [0.5 NumNodes+0.5];
ax.XTick = 1:NumNodes;
ax.YTick = 1:NumNodes;
ax.YTickLabel = Names; 
ax.XTickLabel = Names; 
ax.XTickLabelRotation = 45; 

figure;
scatter(vect(:,1), vect(:,2), '.b');
title('Weight-distance relationship'); xlabel('Distance between ROIs'); ylabel('Weight'); 
end

