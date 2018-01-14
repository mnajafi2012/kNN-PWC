function [training, testing] = draw_random_samples(dataMatrix, N, n)
% draw_10_samples
% dataMatrix: Nx3 matrix; where N = 100 and 3 means 3 columns of two
%             position of xy and one column of labels
% training:   nx3 samples; where n = 10 is asked by the problem as training
%             dataset and 3 is the same as above.

% testing:    (N-n)x3; where n is 10, N is 100 and 3 is the same as above.

% Author: Maryam Najafi
% Date: Sep 28, 2016

% Uniformly-randomly chosen 10 points out of 100; keep as indices I
I = datasample(1:N,n,'Replace',false);
training = dataMatrix(I,:);

testing = [];
for i = 1: size(dataMatrix,1)
    if isempty(I(I == i))
        testing = [testing; dataMatrix(i, :)];
    end 
end

end