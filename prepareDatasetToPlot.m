function training = prepareDatasetToPlot(dataMatrix)
% prepareDatesetToPlot

% Author: Maryam Najafi
% Created Date: Sep 27, 2016

training = cell(2,1);
training{1} = dataMatrix((dataMatrix(:,3) == 1), 1:2); % exclude labels
training{2} = dataMatrix((dataMatrix(:,3) == 2), 1:2); % exclude labels
end