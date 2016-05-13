function dataMatrixNorm = Norm_hampel(dataMatrix)

dataMatrixNorm = zeros(size(dataMatrix));
numFeatures = size(dataMatrix,2);
        for i = 1:numFeatures % cycle through the features
            dataMatrixNorm(:,i) = Normalisation_tanh_hampel(dataMatrix(:,i));
        end

end