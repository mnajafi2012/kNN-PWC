
function plotDecisionBoundary(data, n, classifier)
% function plotDecisionBoundary()
%plotDecisionBoundary
% The following function code is borrowed from this website:
% Ref: http://www.peteryu.ca/tutorials/matlab/visualize_decision_boundaries

% training: Structure of two cells each of which contains a matrix of
%           mx2; where m is a number between 0 and n, and 2 represents xy
%           position of the point on the screen
%           cell1 includes data samples from only class1
%           cell2 includes data samples from only class2


% modified by Maryam Najafi on Sep 27, 2016
% Some explanatory sentences above lines are from the author of the 
% reference link.

% Author: Maryam Najafi
% Created Date: Sep 27, 2016


training = prepareDatasetToPlot(data);

% set up the domain over which you want to visualize the decision
% boundary
xrange = [0 1];
yrange = [0 1];
% step size for how finely you want to visualize the decision boundary.
inc = 0.01;
 
% generate grid coordinates. this will be the basis of the decision
% boundary visualization.
[x, y] = meshgrid(xrange(1):inc:xrange(2), yrange(1):inc:yrange(2));
 
% size of the (x, y) image, which will also be the size of the 
% decision boundary image that is used as the plot background.
image_size = size(x);
 
xy = [x(:) y(:)]; % make (x,y) pairs as a bunch of row vectors.
xy = [reshape(x, image_size(1)*image_size(2),1) reshape(y, image_size(1)*image_size(2),1)];
numxypairs = length(xy); % number of (x,y) pairs
 
% distance measure evaluations for each (x,y) pair.
dist = [];

k = 1; % task 2. a, c)
% p = 1; % task 2. a, b)
% k = 5; % task 2. b)
p = inf; % task 2. c)

kernel_type = 2; % 1: Gaussian, 2: Squared Sinc
spread = 0.02;    % The spread (covariance) of each kernel

global unknown_label
idx = [];
X = [];
% make the background unlabeled matrix X
for i=0:inc:1
    for j = 0:inc:1
        X = [X; [i,j]];
%         idx = [idx; knnClassifier(training, k, vec)];
    end
end

% classify
if (isequal(classifier, 'knn'))
    [Ypred, PCP] = knn_classify(X, data, k, p, unknown_label);
end
if (isequal (classifier, 'pwc'))
    [Ypred, PCP] = pwc_classify( X, data, kernel_type, spread, unknown_label );
end

% reshape the idx (which contains the class label) into an image.
decisionmap = reshape(Ypred, image_size);
figure;
 
%show the image
imagesc(xrange,yrange,decisionmap);
hold on;
set(gca,'ydir','normal');
 
% colormap for the classes:
% class 1 = light red, 2 = light green, 3 = light blue
cmap = [1 0.8 0.8; 0.9 0.9 1];
colormap(cmap);
 
% plot the class training data.
plot(training{1}(:,1),training{1}(:,2), 'ro');
plot(training{2}(:,1),training{2}(:,2), 'bo');
 
% include legend
% legend('Class 1', 'Class 2','Location','NorthOutside', ...
%     'Orientation', 'horizontal');
legend('Class 1', 'Class 2','Location','NorthOutside', ...
    'Orientation', 'horizontal');

if (isequal(classifier, 'knn'))
    title (sprintf('n= %d, k = %d, p = %d', n,k,p));
end
if (isequal (classifier, 'pwc'))
    title (sprintf('n= %d, Kernel Type = %d, spread = %d', n,kernel_type, spread));
end
if (isequal (classifier, 'sqrtSinc'))
    title (sprintf('n= %d, Kernel Type = %d, spread = %d', n,kernel_type, spread));
end

% label the axes.
xlabel('x');
ylabel('y');
hold on;
set(gca,'ydir','normal');

% plot the Optimal decision boundary
drawOptimalDecisionBoundary();
axis equal
xlim ([0 1]);
ylim ([0 1]);

%% Task 4. c
% Plot the posterior class probability in 3D
figure();
global posterior_class_1
global posterior_class_2
base = X(:,2);
tt1 = zeros (image_size(1) , image_size(1));
for i = 1:image_size(1)
    tt1(:, i) = base(i);
end
tt2 = tt1';
plotdata = reshape (posterior_class_1, image_size);
surfc(tt1,tt2, plotdata)
title(sprintf('Posterior Probability Functions for Class 1 \n number of samples = %d', n));
xlabel ('x'); ylabel('y'); zlabel('Posterior Class Probability');
% ('PWC using Gaussian Kernel, spread = 0.1');

% add posterior class 2 to the plot
% hold on
% plotdata2 = reshape (posterior_class_2, image_size);
% surfc(tt1,tt2, plotdata2, 'FaceColor','interp','FaceLighting','gouraud')

end