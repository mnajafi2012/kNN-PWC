function error = cal_error(testing, Ypred)
% cal_knn_error_rate
% testing:  nx3 the rest of randomly chosen data samples. n = 90 by 3 is xy
%           position plus the label
% Ypred:    predicted labels using k-NN for testing data
% error:    The calculated error rate out of knn algorithm

% Author: Maryam Najafi
% Create Date: Sep 28, 2016

n_test = size(testing,1);
error = 0; % the number of misclassified samples
error = length (testing (testing(:,3) ~= Ypred));
error = error / size(testing,1);
end