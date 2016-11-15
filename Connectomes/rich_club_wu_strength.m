function   [Rw] = rich_club_wu_strength(CIJ,varargin)
%RICH_CLUB_WU 	Rich club coefficients curve (weighted undirected graph)
%
%   Rw = rich_club_wu(CIJ,varargin) % rich club curve for weighted graph
%
%   The weighted rich club coefficient, Rw, at level k is the fraction of
%   edge weights that connect nodes of degree k or higher out of the
%   maximum edge weights that such nodes might share.
%
%   Inputs:
%       CIJ:        weighted directed connection matrix
%
%       k-level:    (optional) max level of RC(k).
%                   (by default k-level quals the maximal degree of CIJ)
%                
%   Output:
%       Rw:         rich-club curve
%
%
%   References:     
%       T Opsahl et al. Phys Rev Lett, 2008, 101(16)
%       M van den Heuvel, O Sporns, J Neurosci 2011 31(44)
%
%   Martijn van den Heuvel, University Medical Center Utrecht, 2011

%   Modification History:
%   2011: Original
%   2015: Expanded documentation (Mika Rubinov)


NofNodes = size(CIJ,2);     %#ok<NASGU> %number of nodes
NodeStrength = sum(CIJ); %define strength of each node

%define to which level rc should be computed
if size(varargin,2)==1
    slevel = varargin{1};
elseif isempty(varargin)
    slevel = floor(max(NodeStrength));
else
    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
end

%wrank contains the ranked strength of the network, with strongest nodes on top
strength = sum(CIJ);
srank = sort(strength, 'descend');
smin = min(NodeStrength);

%loop over all possible k-levels
svector = linspace(smin, slevel, 50);
for kk = 1:length(svector)
    
    
    SmallNodes=find(NodeStrength<svector(kk));
    
    if isempty(SmallNodes);
        Rw(kk)=NaN;             %#ok<*AGROW>
        continue
    end
    
    %remove small nodes with NodeStrength<kk
    CutoutCIJ=CIJ;
    CutoutCIJ(SmallNodes,:)=[];
    CutoutCIJ(:, SmallNodes)=[];

    
    %total strength of connections in subset E>r
    Wr = sum(CutoutCIJ(:));
    
    %total number of nodes in subset E>r
    Er = length(sum(CutoutCIJ~=0));
    
    %E>r number of connections with max strength in network
    srank_r = srank(1:1:Er);
    
    %weighted rich-club coefficient
    Rw(kk)=Wr / sum(srank_r); 
end