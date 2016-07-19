function [Rich, Feeder, Local, MeanAll, x] = RichDegreeBFcorrected(Rexp,Adj)
% Rexp - coexpression matrix
% Adj - adjecency matrix

NumNodes = size(Adj, 1);
Degree = degrees_und(Adj);
klevel = max(Degree);
Adjlog = logical(Adj);

Rich = zeros(klevel,1);
Feeder = zeros(klevel,1);
Local = zeros(klevel,1);
MeanAll = zeros(klevel,1);

    for k = 1:klevel
        
        isHub = Degree > k;
        mask = zeros(size(Adj));
        mask(isHub,isHub) = 1; % rich
        mask(isHub,~isHub) = 2; % feeder
        mask(~isHub,isHub) = 2;
        mask(~isHub,~isHub) = 3; % peripheral

        %check no zeros left
        if any(mask(:)==0)
            error('Problem making mask -- zeros remain');
        end
        % Maybe remove diagonal/lower triangle
        %%make a mask for rich, feeder and periferal connection
        
        mask = mask.*Adjlog;
        Rich(k) = mean(Rexp(mask==1));
        Feeder(k) = mean(Rexp(mask==2));
        Local(k) = mean(Rexp(mask==3));
        ToMean = [Rich(k) Feeder(k) Local(k)];
        MeanAll(k) = mean(ToMean);
        
    end

x = linspace(1,klevel,klevel);
figure; plot(x, Rich, x, Feeder, x, Local, x, MeanAll); xlabel('Degree, k'); ylabel('Coexpression Z scored'); legend('Rich', 'Feeder', 'Local', 'MeanALL');
end
