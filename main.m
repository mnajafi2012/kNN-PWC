% Dataset: CIS Noisy Circle-in-the-Square
% The objective is to classify training samples using K-NN and PWC
% algorithms, and compare their misclassification errors with
% Bayes Classifier.

% Author: Maryam Najafi
% Created date: Sep 25, 2016
% Last modified date: Oct 2, 2016

% close all
% clc
% clear all
% global unknown_label
% unknown_label= 0;

% % load ('CIS.mat');
% % n = size(label, 1);
% global D;
% D = 3;
% % figure();
% % scatter (x(label==1), y(label==1), 'r'); hold on; scatter (y(label ==2), x(label ==2), 'b');
% 
% % Generate Uniformly distributed random probabilities
% %   rand_p = rand(n,1);
%    p_flip = 0.5;
% %   label(label(find(rand_p < p_flip)) == 2) = 0;
% %   label(label(find(rand_p < p_flip)) == 1) = 2;
% %   label(label == 0 ) = 1;
% 
% % The noisy data is saved as NCIS.mat file and from now on this file is
% % opened instead of CIS.mat.
% load ('NCIS.mat');
% N = size(label, 1);
% figure(); 
% scatter (x(label==1), y(label==1), 'r'); hold on; scatter (y(label ==2), x(label ==2), 'b');title ('NCIS dataset');

%% Task 1. a) 
% % Prior Probabilities
% prior_1 = size(label(label(:) == 1)) / N;
% prior_2 = 1 - prior_1;
% 
% % Class Conditional Densities
% area_1 = 0.5;
% area_2 = 1 - area_1;
% density_1 = size(label(label==1)) / area_1; density_1 = density_1(1);
% density_2 = size(label(label==2)) / area_2; density_2 = density_2(1);

%% Task 1. d)
% % Bayes Error Rate
% misclassified = 0;
% misclassified = misclassified + size(rand_p(label(rand_p < p_flip) == 1), 1);
% misclassified = misclassified + size(rand_p(label(rand_p > p_flip) == 2), 1);
% % almost half of the entire data has been corrupted.
% error = misclassified / N;

%% Task 1. e)
% % Draw 100 samples out of NCIS and change pflip accordingly. Plot the output.
% CIS = load ('CIS.mat');
% label_100_s = label(1:100);
% x_100_s = x(1:100);
% y_100_s = y(1:100);
% rand_p_100_s = rand_p(1:100);
% for p_flip=0:0.1:0.4
%     figure();
%     label_100_s(label_100_s(find(rand_p_100_s < p_flip)) == 2) = 0;
%     label_100_s(label_100_s(find(rand_p_100_s < p_flip)) == 1) = 2;
%     label_100_s(label_100_s == 0 ) = 1;
%     scatter (x_100_s(label_100_s==1), y_100_s(label_100_s==1), 'r'); hold on; scatter (y_100_s(label_100_s ==2), x_100_s(label_100_s ==2), 'b'); hold on;
%     drawOptimalDecisionBoundary();
% %     plot (x_circle, y_circle);
%     axis equal; xlim([0,1]); ylim([0,1]);
%     title (['p_{flip} = ' num2str(p_flip) ]);
% end
% close all

%% Task 2. a)

% p_flip = 0.1; % asked by the problem
% N = 100; % asked by the problem
% 
% % draw 100 samples randomly (uniform) dataMatrix = draw_samples_NCIS(CIS,
% % p_flip, N); save dataMatrix;
% 
% % load saved dataMatrix from now on
% load dataMatrix
% 
% step = 20;
% n = 10;
% while n <= 100
%     % take a piece of n samples from dataMatrix
%     data = dataMatrix(1:n,:);
%     
%     % within this function I invoked knn_classify function
%     plotDecisionBoundary(data, n, 'knn');
%     
%     if n < 50
%         step = 20;
%     else
%         step = 25;
%     end
%     n = n + step;
% end

%% Task 3. a)
% % Error rate and Box plot
% load dataMatrix
% load ('NCIS.mat');
% N = 1000;
% dataMatrix1000(:,1) = x;
% dataMatrix1000(:,2) = y;
% dataMatrix1000(:,3) = label;
% base = []; b = 1; b_axis = [];% xaxis of the boxplots
% boxes = [];                   % yaxis of the boxplots
% 
% 
% N = size(dataMatrix1000,1);
% k = 5;
% p = 1;
% e_rate = [];
% E_RATE = []; % task 3. b to record e_rate from this task plus task 3. b.
% % as asked in the problem draw 10 samples randomly for 30 times
% for i = 1 : 30
%     n = 10;
%     [training, testing] = draw_random_samples(dataMatrix1000, N, n);
%     
%     % test using k-NN
%     [Ypred, PCP] = knn_classify(testing, training, k, p, unknown_label);
%     
%     % calculate the k-NN error rate
%     e = cal_error(testing, Ypred);
%     t = e / (N - n);
%     e_rate = [e_rate; t];
%     
% end
% 
% % record
% E_RATE = [E_RATE e_rate];
% base = [base; b * ones(size(E_RATE(:,b)))];
% boxes = [boxes; e_rate];
% b_axis = [b_axis, n];
% 
% boxplot(e_rate, 'Labels',{n}); ylim ([0 1]); 
% xlabel ('Number of training samples');
% ylabel('error rate');
% title(sprintf('k-NN error rate for k =%d and p=%d', k,p));

%% Task 3. b)
% % Do as the previous experiment but this time run it w/ different number of
% % training samples
% % keep the previous outcome to plot with the following boxplots
% 
% for a = 1.2:0.2:2.8
%     e_rate = [];
%     b = b + 1;
%     for i = 1 : 30
%         n = ceil(10^a);
%         [training, testing] = draw_random_samples(dataMatrix1000, N, n);
% 
%         % test using k-NN
%         [Ypred, PCP] = knn_classify(testing, training, k, p, unknown_label);
% 
%         % calculate the k-NN error rate
%         e = cal_error(testing, Ypred);
%         t = e / (N - n);
%         e_rate = [e_rate; t];  
%     end
%     
%     % record
%     E_RATE = [E_RATE e_rate];
%     
%     base = [base; b * ones(size(E_RATE(:,b)))];
%     boxes = [boxes; e_rate];
%     b_axis = [b_axis, n];
% end
% 
% boxplot(boxes, base, 'Labels', b_axis);
% ylim ([0 1]); 
% xlabel ('Number of training samples');
% ylabel('error rate');
% title(sprintf('k-NN error rate for k =%d and p=%d', k,p));

%% Task 4. a) Parzen Window Classification with Gaussian Kernels
% clear
% close all
% clc
% 
%  % find it initialized in pwc_classify
% NCIS = load ('dataMatrix.mat');
% dataMatrix1000 = [NCIS.x NCIS.y NCIS.label]; % Nx3 matrix N = 1000
% N = size(dataMatrix1000,1);
% 
% step = 20;
% n = 10;
% while n <= 100
%     
%     % take a piece of n samples from dataMatrix
%     data = dataMatrix1000(1:n,:);
%     
%     % within this function I invoked pwc_classify function
%     plotDecisionBoundary(data, n, 'pwc');
%     % kernel type is defined through the above function using number 1 & 2.
%     
%     if n < 50
%         step = 20;
%     else
%         step = 25;
%     end
%     n = n + step;
% end

%% Task 5. a) Find the champion k-NN classifier

% clc
% clear
% close all
% % 
% % load ('NCIS.mat');
% % 
% % data = [x y label];
% % N = size(data,1);
% % n_train = 30;
% % n_val = 100;
% % n_test = 100;
% % % split the data into three sets of training, validation, testing samples
% % [training, out] = draw_random_samples(data, N , n_train);
% % validation = out(1:n_val, :);
% % testing = out(n_val+1:n_val + n_test, :);

% load ('train_val_test.mat');
% 
% global unknown_label
% unknown_label= 0;
% p = 1;
% 
% % Calculate misclassification error using validation samples
% error = [];
% for k = 1: 15
%     % Use the validation set in the k_NN classification for k=1,...,15
%     [Ypred, PCP] = knn_classify(validation,training,k,p,unknown_label);
%     
%     % calculate the estimated error
%     error = [error; cal_error(validation, Ypred)];
% end
% k = [1:15];
% bar (k,error);
% xlabel('k'); ylabel('Error'); xlim([1 15]); ylim([0 1]);
% title('Estimated Error vs. k');
% e1 = error';
% 
% % Calculate misclassification error using training samples
% figure();
% error2 = [];
% for k = 1: 15
%     % Use the validation set in the k_NN classification for k=1,...,15
%     [Ypred, PCP] = knn_classify(training,training,k,p,unknown_label);
%     
%     % calculate the estimated error
%     error2 = [error2; cal_error(training, Ypred)];
% end
% k = [1:15];
% bar (k,error2);
% xlabel('k'); ylabel('Error'); xlim([1 15]); ylim([0 1]);
% title('Estimated Error vs. k');
% e2 = error2';

%% Task 5. b. Find the optimal kernel spread value of s using PWC of Gaussian kernels
clc
clear
close all

load ('train_val_test.mat');

global unknown_label
unknown_label= 0;
p = 1;
kernel_type = 1;

% Calculate misclassification error using validation samples
error = [];
max_boundary = 0.6;
for s = 0.001: 0.01:max_boundary
    % Use the validation set
    [Ypred, PCP] = pwc_classify( validation, training, kernel_type, s, unknown_label );
    
    % calculate the estimated error
    error = [error; cal_error(validation, Ypred)];
end
s = [0.01:0.01: max_boundary];
bar (s,error);
xlabel('s'); ylabel('Misclassification Error'); xlim([0.01 max_boundary]); ylim([0 1]);
title('Estimated Error vs. s');
max(error);

