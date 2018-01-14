function [ Ypred, PCP ] = pwc_classify( X, Xref, kernel_type, spread, unknown_label )
%pwc_classify
% Estimates the Probability Density based on Parzen Window Classification
%   X: an Nx by D matrix containing Nx D-dimensional test samples in rows
%   Xref: an N by D matrix containing N D-dimensional training samples in
%   rows
%   kernel_type: 1: Gaussian, 2: Squared Sinc
%   spread: the variance of how the kernels are spreaded out
%   unknown_label: an alternative label for undecided labels for samples

%   Ypred: contains all predicted labels

% Author: Maryam Najafi
% Date created: Sep 29, 2016
% Last date modified: Sep 30, 2016

global posterior_class_1
global posterior_class_2
PCP = 0;

% if the testing data has labels remove it temporarily
if size(X,2) == 3
   X = X(:,1:2); 
end

% The number of training samples
N = size(Xref,1);

% split the training set into two classes
X_ref_1 = Xref(find(Xref(:,3)==1),:);
X_ref_2 = Xref(find(Xref(:,3)==2),:);
N_1 = size(X_ref_1, 1);
N_2 = size(X_ref_2, 1);

Ypred = [];
posterior_class_1 = [];
SIGMA_1 = 0; % The Parzen Probability Density Function Estimate for test samples being in class 1
SIGMA_2 = 0; % The Parzen Probability Density Function Estimate for test samples being in class 2
if (kernel_type == 1) % 1: Gaussian
    for i = 1:size (X, 1)
    candidate = X(i, :);
    
    for j = 1:N_1
        ref = X_ref_1(j, 1:2);
        SIGMA_1 = SIGMA_1 + (1/ ((sqrt(2 * pi) * spread) )) * exp(- (norm((ref - candidate), 2)^2/( 2 * spread^2))); 
    end
    % get the average of Parzen PDF Estimate at point candidate
    SIGMA_1 = SIGMA_1 / N_1;
    
    %% record posterior probabilities for class 1 to plot in 3D in task 4.c.
    posterior_class_1 = [posterior_class_1; SIGMA_1];
    %%
    for j = 1:N_2
        ref = X_ref_2(j, 1:2);
         SIGMA_2 = SIGMA_2 + (1/ ((sqrt(2 * pi) * spread) )) * exp(- (norm((ref -  candidate), 2)^2/ (2 * spread^2))); 
    end
    % get the average of Parzen PDF Estimate at point candidate
    SIGMA_2 = SIGMA_2 / N_2;
    
    %%
     posterior_class_2 = [posterior_class_2; SIGMA_2];
    %%
    
    
    % now check which class the candidate belongs to
    % the one with the higher value for sigma is the ultimate label.
    if (SIGMA_1 < SIGMA_2)
       % means candidate belongs to class 2
       Ypred = [Ypred; 2];
    elseif (SIGMA_1 > SIGMA_2)
       % means candidate belongs to class 1 
       Ypred = [Ypred; 1];
    else
        % means we are on the boundary (make the label unknown-label)
        Ypred = [Ypred; unknown_label];
    end
    
    end

elseif (kernel_type == 2)
SIGMA_1 = 0; % The Parzen Probability Density Function Estimate for test samples being in class 1
SIGMA_2 = 0; % The Parzen Probability Density Function Estimate for test samples being in class 2

    for i = 1:size (X, 1)
    candidate = X(i, :);
    
    for j = 1:N_1
        ref = X_ref_1(j, 1:2);
        u = norm (ref - candidate)/spread;
        if (u == 0)
            SIGMA_1 = SIGMA_1 + (1/pi);
        else 
            SIGMA_1 = SIGMA_1 + (sin(u) / u) ^2;
        end
        
    end
    % get the average of Parzen PDF Estimate at point candidate
    SIGMA_1 = SIGMA_1 / N_1;
    
    for j = 1:N_2
        ref = X_ref_2(j, 1:2);
        u = norm (ref - candidate)/spread;
        if (u == 0)
            SIGMA_2 = SIGMA_2 + (1/pi);
        else 
            SIGMA_2 = SIGMA_2 + (sin(u) / u) ^2;
        end
    end
    % get the average of Parzen PDF Estimate at point candidate
    SIGMA_2 = SIGMA_2 / N_2;
    
    % now check which class the candidate belongs to
    % the one with the higher value for sigma is the ultimate label.
    if (SIGMA_1 < SIGMA_2)
       % means candidate belongs to class 2
       Ypred = [Ypred; 2];
    elseif (SIGMA_1 > SIGMA_2)
       % means candidate belongs to class 1 
       Ypred = [Ypred; 1];
    else
        % means we are on the boundary (make the label unknown-label)
        Ypred = [Ypred; unknown_label];
    end
    
    end

    
end
end

