%Separate distance correction for left hemisphere, right hemisphere and for
%cortex-cortex,cortex-subcortex, subcortex-subcortec

%% not finished script, waiting for final parcallations
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

load(sprintf('%d_DistThresh2_%s_combined_WholeBrainNormalised.mat', NumNodes, Parcellation{1}));

[sROIs, ind] = sort(ROIs);
ExpressionNormSORTED = ExpressionNorm(ind, ind);
MRIvoxCoordinatesSORTED = MRIvoxCoordinates(ind, ind);

BrainPart = {'Cortex', 'Subcortex','Mix'};     
%'LCC', 'LCS', 'LSS', 'RCC', 'RCS', 'RSS','RLCC', 'RLCS1', 'RLCS2', 'RLSS', 'RLCC', 'RLCS1', 'RLCS2'
for part = BrainPart
    
    switch part{1}

        case 'Cortex'
            Fit = {'exp'};

            %left cortex
           [ROIsLC] = find(sROIs<=Lcortex);
            ExpLCC = ExpressionNormSORTED(ROIsLC, ROIsLC); 
            ExpLCC(logical(eye(size(ExpLCC)))) = 0;
            
            
            DistLCC = MRIvoxCoordinatesSORTED(ROIsLC, ROIsLC);
            DistLCC(DistLCC==0) = 0.1;
            DistLCC(logical(eye(size(DistLCC)))) = 0;
            
            % right cortex
            [ROIsRCC] = find(sROIs>Lsubcortex & sROIs<=Rcortex);
            ExpRCC = ExpressionNormSORTED(ROIsRCC, ROIsRCC); 
            ExpRCC(logical(eye(size(ExpRCC)))) = 0;
 

            DistRCC = MRIvoxCoordinatesSORTED(ROIsRCC, ROIsRCC);
            DistRCC(DistRCC==0) = 0.1;
            DistRCC(logical(eye(size(DistRCC)))) = 0;
            
            
            % left - right 
            [ROIsRLCC1] = find(sROIs<=Lcortex);
            [ROIsRLCC2] = find(sROIs>Lsubcortex & sROIs<=Rcortex);
            ExpRLCC = ExpressionNormSORTED(ROIsRLCC2, ROIsRLCC1); 
            DistRLCC = MRIvoxCoordinatesSORTED(ROIsRLCC2, ROIsRLCC1); 
            
            
            M1Exp = vertcat(ExpLCC, ExpRLCC);
            M2Exp = vertcat(ExpRLCC', ExpRCC);
            MExp = horzcat(M1Exp, M2Exp); 
            
            M1Dist = vertcat(DistLCC, DistRLCC);
            M2Dist = vertcat(DistRLCC', DistRCC);
            MDist = horzcat(M1Dist, M2Dist);
       
            Rvect = MExp(:); 
            indR = logical(Rvect(:))+0;
            Rvect = nonzeros(MExp(:)); 
            indR(indR>0) = Rvect; 
            
            
            
            
            Dvect = MDist(:); 
            indD = logical(Dvect(:))+0;
            indF = logical(Dvect(:))+0;
            Dvect = nonzeros(MDist(:)); 
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
            NumSamples = size(ExpLCC,1);
            ExpLCC = reshape(Residuals,[NumSamples, NumSamples]);

            BF_PlotQuantiles(Dvect,(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
            
        
        case 'LSS'

            
            Fit = {'exp'};

            [ROIsLS] = find(sROIs>Lcortex & sROIs<=Lsubcortex);
            ExpLSS = ExpressionNormSORTED(ROIsLS, ROIsLS); 
            ExpLSS(logical(eye(size(ExpLSS)))) = 0;
            
            Rvect = ExpLSS(:); 
            indR = logical(Rvect(:))+0;
            
            Rvect = nonzeros(ExpLSS(:)); 
            
            indR(indR>0) = Rvect; 

            DistLSS = MRIvoxCoordinatesSORTED(ROIsLS, ROIsLS);
            DistLSS(DistLSS==0) = 0.1;
            DistLSS(logical(eye(size(DistLSS)))) = 0;
            
            Dvect = DistLSS(:); 
            indD = logical(Dvect(:))+0;
            indF = logical(Dvect(:))+0;
            
            Dvect = nonzeros(DistLSS(:)); 
            
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
            NumSamples = size(ExpLSS,1);
            ExpLSS = reshape(Residuals,[NumSamples, NumSamples]);

            BF_PlotQuantiles(Dvect,(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
            
            
%             Fit = {'exp'};
% 
%             [ROIsLS] = find(sROIs>Lcortex & sROIs<=Lsubcortex);
%             ExpLSS = ExpressionNormSORTED(ROIsLS, ROIsLS); 
%             ExpLSS(logical(eye(size(ExpLSS)))) = 0;
%             DistLSS = MRIvoxCoordinatesSORTED(ROIsLS, ROIsLS); 
%             DistLSS(logical(eye(size(DistLSS)))) = 0;
%             
%             
% 
%             Rvect = ExpLSS(:); 
%             Dvect = DistLSS(:);
%             
%             
%             BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
%             [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
%             FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
%             plot(c); 
% 
%             Residuals = Rvect - FitCurve;
%             NumSamples = size(ExpLSS,1);
%             ExpLSS = reshape(Residuals,[NumSamples, NumSamples]);
% 
%             BF_PlotQuantiles(Dvect,nonzeros(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
        case 'LCS'
        %% cortex - subcortex

            Fit = {'exp'};

            [ROIsLCS1] = find(sROIs>Lcortex & sROIs<=Lsubcortex);
            [ROIsLCS2] = find(sROIs<=Lcortex);
            ExpLCS = ExpressionNormSORTED(ROIsLCS1, ROIsLCS2); 
            %ExpLCS(logical(eye(size(ExpLCS)))) = 0;
            DistLCS = MRIvoxCoordinatesSORTED(ROIsLCS1, ROIsLCS2); 
            %DistLCS(logical(eye(size(DistLCS)))) = 0;

            Rvect = ExpLCS(:); 
            Dvect = DistLCS(:);
            
            
            BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
            [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
            FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
            plot(c); 

            Residuals = Rvect - FitCurve;
            NumSamples1 = size(ExpLCS,1);
            NumSamples2 = size(ExpLCS,2);
            ExpLCS = reshape(Residuals,[NumSamples1, NumSamples2]);

            BF_PlotQuantiles(Dvect,nonzeros(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
            
            case 'RCC'
        %% right cortex - cortex
            
             Fit = {'exp'};

            [ROIsRCC] = find(sROIs>Lsubcortex & sROIs<=Rcortex);
            ExpRCC = ExpressionNormSORTED(ROIsRCC, ROIsRCC); 
            ExpRCC(logical(eye(size(ExpRCC)))) = 0;
            
            Rvect = ExpRCC(:); 
            indR = logical(Rvect(:))+0;
            
            Rvect = nonzeros(ExpRCC(:)); 
            
            indR(indR>0) = Rvect; 

            DistRCC = MRIvoxCoordinatesSORTED(ROIsRCC, ROIsRCC);
            DistRCC(DistRCC==0) = 0.1;
            DistRCC(logical(eye(size(DistRCC)))) = 0;
            
            Dvect = DistRCC(:); 
            indD = logical(Dvect(:))+0;
            indF = logical(Dvect(:))+0;
            
            Dvect = nonzeros(DistRCC(:)); 
            
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
            NumSamples = size(ExpRCC,1);
            ExpRCC = reshape(Residuals,[NumSamples, NumSamples]);

            BF_PlotQuantiles(Dvect,(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));

%             [ROIsRCC] = find(sROIs>Lsubcortex & sROIs<=Rcortex);
%             ExpRCC = ExpressionNormSORTED(ROIsRCC, ROIsRCC); 
%             %ExpRCC(logical(eye(size(ExpRCC)))) = NaN;
%             DistRCC = MRIvoxCoordinatesSORTED(ROIsRCC, ROIsRCC);
%             %DistRCC(logical(eye(size(DistRCC)))) = NaN;
% 
%             Rvect = ExpRCC(:); 
%             Dvect = DistRCC(:);
%             
%             BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
%             [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
%             FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
%             plot(c); 
% 
%             Residuals = Rvect - FitCurve;
%             NumSamples = size(ExpRCC,1);
%             ExpRCC = reshape(Residuals,[NumSamples, NumSamples]);
% 
%             BF_PlotQuantiles(Dvect,nonzeros(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
            
            
            case 'RSS'
        %% right subcortex - subcortex

            Fit = {'exp'};
            
            [ROIsRSS] = find(sROIs>Rcortex & sROIs<=NumNodes);
            ExpRSS = ExpressionNormSORTED(ROIsRSS, ROIsRSS); 
            ExpRSS(logical(eye(size(ExpRSS)))) = 0;
            
            Rvect = ExpRSS(:); 
            indR = logical(Rvect(:))+0;
            
            Rvect = nonzeros(ExpRSS(:)); 
            
            indR(indR>0) = Rvect; 

            DistRSS = MRIvoxCoordinatesSORTED(ROIsRSS, ROIsRSS);
            DistRSS(DistRSS==0) = 0.1;
            DistRSS(logical(eye(size(DistRSS)))) = 0;
            
            Dvect = DistRSS(:); 
            indD = logical(Dvect(:))+0;
            indF = logical(Dvect(:))+0;
            
            Dvect = nonzeros(DistRSS(:)); 
            
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
            NumSamples = size(ExpRSS,1);
            ExpRSS = reshape(Residuals,[NumSamples, NumSamples]);

            BF_PlotQuantiles(Dvect,(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));

%             [ROIsRSS] = find(sROIs>Rcortex & sROIs<=NumNodes);
%             ExpRSS = ExpressionNormSORTED(ROIsRSS, ROIsRSS); 
%             %ExpRSS(logical(eye(size(ExpRSS)))) = NaN;
%             DistRSS = MRIvoxCoordinatesSORTED(ROIsRSS, ROIsRSS);
%             %DistRSS(logical(eye(size(DistRSS)))) = NaN;
% 
%             Rvect = ExpRSS(:); 
%             Dvect = DistRSS(:);
%             
%             BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
%             [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
%             FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
%             plot(c); 
% 
%             Residuals = Rvect - FitCurve;
%             NumSamples = size(ExpRSS,1);
%             ExpRSS = reshape(Residuals,[NumSamples, NumSamples]);
% 
%             BF_PlotQuantiles(Dvect,nonzeros(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
%             
            case 'RCS'
        %% right subcortex - subcortex

            Fit = {'exp'};

            [ROIsRCS1] = find(sROIs>Lsubcortex & sROIs<=Rcortex);
            [ROIsRCC2] = find(sROIs>Rcortex & sROIs<=NumNodes);
            ExpRCS = ExpressionNormSORTED(ROIsRCS1, ROIsRCC2); 
            %ExpRCS(logical(eye(size(ExpRCS)))) = 0;
            DistRCS = MRIvoxCoordinatesSORTED(ROIsRCS1, ROIsRCC2);
            %DistRCS(logical(eye(size(DistRCS)))) = 0;

            Rvect = ExpRCS(:); 
            Dvect = DistRCS(:);
            
            BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
            [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
            FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
            plot(c); 

            Residuals = Rvect - FitCurve;
            NumSamples1 = size(ExpRCS,1);
            NumSamples2 = size(ExpRCS,2);
            ExpRSS = reshape(Residuals,[NumSamples1, NumSamples2]);

            BF_PlotQuantiles(Dvect,nonzeros(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
            
            
            case 'RLCC'
        %% right left cortex - cortex

            Fit = {'exp'};

            [ROIsRLCC1] = find(sROIs<=Lcortex);
            [ROIsRLCC2] = find(sROIs>Lsubcortex & sROIs<=Rcortex);
            ExpRLCC = ExpressionNormSORTED(ROIsRLCC2, ROIsRLCC1); 
            %ExpRLCC(logical(eye(size(ExpRLCC)))) = 0;
            DistRLCC = MRIvoxCoordinatesSORTED(ROIsRLCC2, ROIsRLCC1); 
            %DistRLCC(logical(eye(size(DistRLCC)))) = 0;

            Rvect = ExpRLCC(:); 
            Dvect = DistRLCC(:);
            
             BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
            [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
            FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
            plot(c); 

            Residuals = Rvect - FitCurve;
            NumSamples1 = size(ExpRLCC,1);
            NumSamples2 = size(ExpRLCC,2);
           ExpRLCC = reshape(Residuals,[NumSamples1, NumSamples2]);

            BF_PlotQuantiles(Dvect,nonzeros(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
            
            
            
            case 'RLSS'
        %% right subcortex - subcortex

            Fit = {'exp'};

            [ROIsRLSS1] = find(sROIs>Lcortex & sROIs<=Lsubcortex);
            [ROIsRLSS2] = find(sROIs>Rcortex & sROIs<=NumNodes);
            ExpRLSS = ExpressionNormSORTED(ROIsRLSS2, ROIsRLSS1); 
            %ExpRLSS(logical(eye(size(ExpRLSS)))) = 0;
            DistRLSS = MRIvoxCoordinatesSORTED(ROIsRLSS2, ROIsRLSS1);
            %DistRLSS(logical(eye(size(DistRLSS)))) = 0;

            Rvect = ExpRLSS(:); 
            Dvect = DistRLSS(:);
            
            BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
            [f_handle,Stats,c] = GiveMeFit(Dvect,Rvect,Fit{1});
            FitCurve = feval('fitCurve',Dvect, Fit, f_handle,Stats,c);
            plot(c); 

            Residuals = Rvect - FitCurve;
            NumSamples1 = size(ExpRLSS,1);
            NumSamples2 = size(ExpRLSS,2);

           ExpRLSS = reshape(Residuals,[NumSamples1, NumSamples2]);

            BF_PlotQuantiles(Dvect,nonzeros(Residuals(:)),20,0,1); title(sprintf('%s Coexpresion vs distance residuals', part{1}));
            
            
            
            
            case 'RLCS1'
        %% right subcortex - subcortex

            Fit = {'exp'};

            [ROIsRLCS11] = find(sROIs>Lsubcortex & sROIs<=Rcortex);
            [ROIsRLCS22] = find(sROIs>Lcortex & sROIs<=Lsubcortex);
            ExpRLCS1 = ExpressionNormSORTED(ROIsRLCS11, ROIsRLCS22); 
            %ExpRLCS1(logical(eye(size(ExpRLCS1)))) = 0;
            DistRLCS1 = MRIvoxCoordinatesSORTED(ROIsRLCS11, ROIsRLCS22); 
            %DistRLCS1(logical(eye(size(DistRLCS1)))) = 0;

            Rvect = ExpRLCS1(:); 
            Dvect = DistRLCS1(:);
            
            BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
            % do not fit anything to RL cortex-subcortex samples
            
            
             case 'RLCS2'
        %% right subcortex - subcortex

            Fit = {'exp'};

            [ROIsRLCS21] = find(sROIs>Rcortex & sROIs<=NumNodes);
            [ROIsRLCS22] = find(sROIs<=Lcortex);
            ExpRLCS2 = ExpressionNormSORTED(ROIsRLCS21, ROIsRLCS22); 
            %ExpRLCS2(logical(eye(size(ExpRLCS2)))) = 0;
            DistRLCS2 = MRIvoxCoordinatesSORTED(ROIsRLCS21, ROIsRLCS22); 
            %DistRLCS2(logical(eye(size(DistRLCS2)))) = 0;

            Rvect = ExpRLCS2(:); 
            Dvect = DistRLCS2(:);
            
            BF_PlotQuantiles(Dvect,Rvect,200,0,1);title(sprintf('%s Coexpresion vs distance fit', part{1}));
            % do not fit anything to RL cortex-subcortex samples
            
    end
end


M1Exp = horzcat( ExpLCC, ExpLCS');
M2Exp = horzcat (ExpLCS, ExpLSS);
ML = cat(1, M1Exp, M2Exp);

M3 = cat(2, ExpRCC, ExpRCS);
M4 = cat(2, ExpRCS', ExpRSS);
MR = cat(1, M3, M4);
% 
% 
M5 = horzcat(ExpRLCC, ExpRLCS1); 
M6 = horzcat(ExpRLCS2, ExpRLSS);
MRL = vertcat(M5, M6);
MRL2 = MRL'; 

% stack parts together

MU = horzcat(ML, MRL2);
MD = horzcat(MRL, MR);
MExp = vertcat(MU, MD);

figure; imagesc(MExp); caxis([-1,1]);title('Sample-sample coexpression corrected');
colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);


W = unique(sROIs);
ParcelCoexpression = zeros(length(W),length(W));
for i=1:length(W)
    for j=1:length(W)

        A = find(sROIs == W(i));
        B = find(sROIs == W(j));
        P = MExp(A, B);

        ParcelCoexpression(i,j) = mean(mean(P));

    end
end

figure; imagesc(ParcelCoexpression); caxis([-1,1]);title('Sample-sample coexpression corrected');
colormap([flipud(BF_getcmap('blues',9));BF_getcmap('reds',9)]);



