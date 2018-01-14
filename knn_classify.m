function [ Ypred, PCP ] = knn_classify( X, Xref, k, p, unknown_label)
%knn_classify classifies each new data entry and assign it to a class.
%According to the problem followings are definitions of function variables.
%   X: an Nx by D matrix containing Nx D-dimensional test samples in rows
%   Xref: an N by D matrix containing N D-dimensional training samples in
%   rows
%   k: the number of nearest-neighbors to be used
%   p: specifies the Lp norm to be used for measuring distances
%   unknown_label: an integer that substitute an undecided label due to
%   voting ties

%   Ypred: contains all predicted labels

% Author: Maryam Najafi
% Date created: Sep 25, 2016
% Last date modified: Sep 27, 2016


% if the testing data has labels remove it temporarily
if size(X,2) == 3
   X = X(:,1:2); 
end



Ypred = [];
PCP = 0;
for i = 1: size (X,1)
    candidate = X(i, :);
    stack = []; % nX2 matrix comprises with column 1: dist 
                                          % column 2: the training sample label
    for j = 1: size(Xref,1)
        ref = Xref(j, 1:2);
        
        dist = norm((candidate - ref), p); % column 1
%         dist = sqrt(sum((candidate - ref) .^ 2));
        label = Xref(j, 3);          % column 2
        
        s = [dist, label];
        stack = [stack; s];
    end
    % sort for the n (e.g. n = 10) distances along with the labels.
    [y,I] = sort(stack(:,1));
    stack = stack(I,:);
    
    % find the first k nearest distances
    labels = stack(1:k, 2);
    num_of_label_1 = length(labels(labels == 1));
    num_of_label_2 = length(labels(labels == 2));
    
    % compute the majority
    if num_of_label_1 < num_of_label_2
        Ypred = [Ypred; 2];
    else if num_of_label_1 > num_of_label_2
            Ypred = [Ypred; 1];
        else
            Ypred = [Ypred; unknown_label];
        end
    end
end

end

