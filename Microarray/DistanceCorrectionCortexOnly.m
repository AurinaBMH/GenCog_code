close all; clear all; 
cd ('/Users/Aurina/Documents/Genetics_connectome/Gen_Cog/Data/Microarray/S01-S06_combined/');

%% left side
%% cortex-cortex

Parcellation = {'cust250'}; 

if strcmp(Parcellation, 'aparcaseg')
                Lcortex = 34;
                Lsubcortex = 41;
                Rcortex = 75;
                NumNodes = 82;
    elseif strcmp(Parcellation, 'cust100')
                NumNodes = 220;
                Lcortex = 100;
                Lsubcortex = 110;
                Rcortex = 210;
    elseif strcmp(Parcellation, 'cust250')
                NumNodes = 530;
                Lcortex = 250;
                Lsubcortex = 265;
                Rcortex = 515;

end

load(sprintf('%d_DistThresh2_%s_combined_CortexNormalised.mat', NumNodes, Parcellation{1}));

[sROIs, ind] = sort(ROIs);
ExpressionNormSORTED = ExpressionNorm(ind, ind);
MRIvoxCoordinatesSORTED = MRIvoxCoordinates(ind, ind);

part = {'Left + Right cortex'};

Fit = {'exp'};


            ExpCortex = ExpressionNormSORTED; 
            ExpCortex(logical(eye(size(ExpCortex)))) = 0;
            
            Rvect = ExpCortex(:); 
            indR = logical(Rvect(:))+0;
            
            Rvect = nonzeros(ExpCortex(:)); 
            
            indR(indR>0) = Rvect; 

            DistLCC = MRIvoxCoordinatesSORTED;
            DistLCC(DistLCC==0) = 0.1;
            DistLCC(logical(eye(size(DistLCC)))) = 0;
            
            Dvect = DistLCC(:); 
            indD = logical(Dvect(:))+0;
            indF = logical(Dvect(:))+0;
            
            Dvect = nonzeros(DistLCC(:)); 
            
            indD(indD>0) = Dvect; 

            
            BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
            [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
            FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
            plot(c); 

            indF(indF>0) = FitCurve;
            Rvect = indR;
            Dvect = indD;
            FitCurve = indF;
            
            Residuals = Rvect - FitCurve;
            NumSamples = size(ExpCortex,1);
            ExpCortex = reshape(Residuals,[NumSamples, NumSamples]);

            BF_PlotQuantiles(Dvect,(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));


% average samples according to ROI
W = unique(sROIs);
ParcelCoexpression = zeros(length(W),length(W));
for i=1:length(W)
    for j=1:length(W)

        A = find(sROIs == W(i));
        B = find(sROIs == W(j));
        P = ExpCortex(A, B);

        ParcelCoexpression(i,j) = mean(mean(P));

    end
end

figure; imagesc(ParcelCoexpression); caxis([-1,1]);title('Sample-sample coexpression corrected');
colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);
