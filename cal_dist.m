function [training] = cal_dist (dist, data, n)
% cal_dist
%   dist: a square matrix nxn where n is given by the problem
%         e.g. n = 10; dist is 10x10
%   dataMatrix: NCIS dataset 100x3

% training = dataMatrix(1:n,:);

% Author: Maryam Najafi
% Created Date: Sep 26, 2016

for i = 1:n
    candidate_vec = data(i, 1:2);
    min = inf;
    for j = 1:n
        if i ~= j
            examin_vec = data(j, 1:2);
            dist(i,j) = sqrt(sum((candidate_vec - examin_vec).^2));
            if (dist(i,j) < min)
               min = dist(i,j);
               training(i,3) = data(j, 3);
               training(i,1:2) = candidate_vec;
            end
%             dataMatrix(i,3) = dataMatrix(min(dist(i,:)),3);
        end
    end
end


end