function [ samples ] = draw_samples_NCIS( dataset, p_flip, N )
%draw_samples_NCIS 
%   Draws 100 samples from CIS and flips labels according to the value of
%   p_flip which is 0.1 in our example.

%   dataset: CIS dataset
%   p_flip: 0.1
%   samples: an NxD data matrix where N = 100 is the number of samples 
%   and D = 3 containing x,y and third column are positions and label.
%   N = 100 the number of samples

% Author: Maryam Najafi
% Created Date: sep 26, 2016

D = 3;

dataset.label(dataset.label(find(dataset.rand_p < p_flip)) == 2) = 0;
dataset.label(dataset.label(find(dataset.rand_p < p_flip)) == 1) = 2;
dataset.label(dataset.label == 0 ) = 1;

samples = zeros(N, D);
dataset.label(1:N);

% Uniformly-randomly chosen 100 points out of 1000
r = datasample(1:1000,N,'Replace',false);

samples(:,1) = dataset.x(r);
samples(:,2) = dataset.y(r);
samples(:,3) = dataset.label(r);

end

