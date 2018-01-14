function drawOptimalDecisionBoundary()
% drawOptimalDecisionBoundary

% Author: Maryam Najafi
% Created Date: Sep 28, 2016

% draw optimal decision boundary
x_origin = 0.5; y_origin = 0.5; r = 0.399;
theta = linspace (0, 2*pi);
x_circle = r * cos (theta) + x_origin; y_circle = r * sin(theta) + y_origin;
plot (x_circle, y_circle);
end