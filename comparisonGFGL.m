%% Comparison with GFGL
% 
%   This file is obtained adapting normalExample.m by
%   Alex Gibberd - UCL Department of Statistical Science
%   from the supplementary material of "Gibberd and Nelson (2017), 
%   Regularized estimation of piecewise constant Gaussian graphical models
%   the group-fused graphical lasso."
%   Moreover the file use the fuction clustCoeff_modified, which we adapted
%   from 
%   github.com/ivanbrugere/matlab-networks-toolbox
%   to compute Global clustering coefficients.

%   TO USE THIS CODE, you need to install the MATLAB code contained in the
%   supplementary of Gibberd and Nelson (2017). 

%   WARNING: 
%   Gibberd and Nelson define the change point (CP_GB) as a point
%   after which the covariance changes. We define the change point as 
%   the first time point when the change is observed. 
%   Thus, for us, CP = CP_GB + 1
%   Plots of change points here are accordingly to our definition, however
%   the vector of change points in the code, i.e. cpG, is encoded
%   accordlingly to the definition in Gibberd and Nelson (2017)


%% Set Default parameters
warning('off'); % Overwrite warnings when viewing via publish
lambda1G=0.25;
lambda2G=20;
%% Simulated Data - scenario 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%code for 1 replica
y = csvread("C:\Users\beatr\Documents\SimData_scenario3\Ys3r1.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.
display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);
%%TRUe GRAPH STRUCTURE 
AdTrue1 = zeros(10); 
for i=1:9
    AdTrue1(i,i+1) = 1;
    AdTrue1(i+1,i) = 1;
    Ad1(i,i) = 0;
    Ad2(i,i) = 0;
end
AdTrue2 = zeros(10); 
AdTrue2(2,9) = 1;  AdTrue2(9,2) = 1;
AdTrue2(8,9) = 1;  AdTrue2(9,8) = 1;
AdTrue2(8,3) = 1;  AdTrue2(3,8) = 1;
AdTrue2(3,4) = 1;  AdTrue2(4,3) = 1;
AdTrue2(4,7) = 1;  AdTrue2(7,4) = 1;
AdTrue2(7,6) = 1;  AdTrue2(6,7) = 1;
AdTrue2(6,5) = 1;  AdTrue2(5,6) = 1;
AdTrue2(6,8) = 1;  AdTrue2(8,6) = 1;
AdTrue2(5,10) = 1;  AdTrue2(10,5) = 1;

TP = 0; FN = 0; TN = 0; FP = 0;
lambda1G=0.5;
lambda2G=60;
%code for 1 replica

y = csvread("C:\Users\beatr\Documents\SimData_scenario3\Ys3r1.csv")
display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

Ad1 = ZG(:,:,35)~=0;
Ad2 = ZG(:,:,135)~=0;


 TP = TP + sum(sum(Ad1 .* AdTrue1))/2 + sum(sum(Ad2 .* AdTrue2))/2;
 FN = FN + sum(sum((1-Ad1) .* AdTrue1))/2 + sum(sum((1-Ad2) .* AdTrue2))/2;
 TN = TN + (sum(sum((1-Ad1) .* (1-AdTrue1)))-10 + sum(sum((1-Ad2) .* (1-AdTrue2)))-10)/2;
 FP = FP + sum(sum((Ad1) .* (1-AdTrue1)))/2 + sum(sum((Ad2) .* (1-AdTrue2)))/2;

 FPR = FP / (FP + TN)
 TPR = TP / (TP + FN)

figure(1)
set(gcf, 'Color', [1,1,1]);
subplot(1,2,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,2,2);
plotGraph(squeeze(abs(ZG(:,:,69))),9,0);
exportgraphics(gcf,'G_scenario1_R1_hyper2.png','Resolution',300)

%% Simulated Data - scenario 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = csvread("C:\Users\beatr\Documents\Y_scenario4.csv")

%hyperpar 1
lambda1G=0.2;
lambda2G=60;
display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_scenario4_0.5.png','Resolution',300)

%% Graphical estimation

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,4,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,4,2);
plotGraph(squeeze(abs(ZG(:,:,(60)))),9,0);
subplot(1,4,3);
plotGraph(squeeze(abs(ZG(:,:,(141)))),9,0);
subplot(1,4,4);
plotGraph(squeeze(abs(ZG(:,:,(150)))),9,0);
exportgraphics(gcf,'G_scenario4_0.5.png','Resolution',300)

AdTrue1 = zeros(20); 
AdTrue1(4,16) = 1;  AdTrue1(16,4) = 1;
AdTrue1(16,6) = 1;  AdTrue1(6,16) = 1;
AdTrue1(6,7) = 1;  AdTrue1(7,6) = 1;
AdTrue1(7,3) = 1;  AdTrue1(3,7) = 1;
AdTrue1(3,18) = 1;  AdTrue1(18,3) = 1;
AdTrue1(18,20) = 1;  AdTrue1(20,18) = 1;
AdTrue1(2,14) = 1;  AdTrue1(14,2) = 1;
AdTrue1(5,1) = 1;  AdTrue1(5,1) = 1;
AdTrue1(17,15) = 1;  AdTrue1(15,17) = 1;
AdTrue1(1,17) = 1;  AdTrue1(17,1) = 1;
AdTrue1(15,13) = 1;  AdTrue1(13,15) = 1;

AdTrue2 = zeros(20); 
AdTrue2(20,18) = 1;  AdTrue2(18,20) = 1;
AdTrue2(18,3) = 1;  AdTrue2(3,18) = 1;
AdTrue2(3,19) = 1;  AdTrue2(19,3) = 1;
AdTrue2(3,7) = 1;  AdTrue2(7,3) = 1;
AdTrue2(7,12) = 1;  AdTrue2(12,7) = 1;
AdTrue2(7,6) = 1;  AdTrue2(6,7) = 1;
AdTrue2(6,16) = 1;  AdTrue2(16,6) = 1;
AdTrue2(16,4) = 1;  AdTrue2(4,16) = 1;
AdTrue2(4,14) = 1;  AdTrue2(14,4) = 1;
AdTrue2(14,2) = 1;  AdTrue2(2,14) = 1;
AdTrue2(2,13) = 1;  AdTrue2(13,2) = 1;
AdTrue2(2,15) = 1;  AdTrue2(15,2) = 1;
AdTrue2(13,15) = 1;  AdTrue2(15,13) = 1;
AdTrue2(15,17) = 1;  AdTrue2(17,15) = 1;
AdTrue2(17,1) = 1;  AdTrue2(17,1) = 1;
AdTrue2(1,5) = 1;  AdTrue2(5,1) = 1;
AdTrue2(5,9) = 1;  AdTrue2(9,5) = 1;

AdTrue3 = zeros(20); 
AdTrue3(9,5) = 1;  AdTrue3(5,9) = 1;
AdTrue3(5,1) = 1;  AdTrue3(1,5) = 1;
AdTrue3(1,17) = 1;  AdTrue3(17,1) = 1;
AdTrue3(17,15) = 1;  AdTrue3(15,17) = 1;
AdTrue3(15,13) = 1;  AdTrue3(13,15) = 1;
AdTrue3(15,14) = 1;  AdTrue3(14,15) = 1;
AdTrue3(15,2) = 1;  AdTrue3(2,15) = 1;
AdTrue3(4,17) = 1;  AdTrue3(17,4) = 1;
AdTrue3(4,14) = 1;  AdTrue3(14,4) = 1;
AdTrue3(4,11) = 1;  AdTrue3(11,4) = 1;
AdTrue3(4,16) = 1;  AdTrue3(16,4) = 1;
AdTrue3(14,2) = 1;  AdTrue3(2,14) = 1;
AdTrue3(13,2) = 1;  AdTrue3(2,13) = 1;
AdTrue3(2,19) = 1;  AdTrue3(19,2) = 1;
AdTrue3(19,3) = 1;  AdTrue3(3,19) = 1;
AdTrue3(3,7) = 1;  AdTrue3(7,3) = 1;
AdTrue3(3,18) = 1;  AdTrue3(18,3) = 1;
AdTrue3(3,12) = 1;  AdTrue3(12,3) = 1;
AdTrue3(12,7) = 1;  AdTrue3(7,12) = 1;
AdTrue3(18,20) = 1;  AdTrue3(20,18) = 1;
AdTrue3(7,6) = 1;  AdTrue3(6,7) = 1;
AdTrue3(6,16) = 1;  AdTrue3(16,6) = 1;
TP = 0; FN = 0; TN = 0; FP = 0;

AdTrue4 = zeros(20); 
AdTrue4(12,3) = 1;  AdTrue4(3,12) = 1;
AdTrue4(7,6) = 1;  AdTrue4(6,7) = 1;
AdTrue4(7,8) = 1;  AdTrue4(8,7) = 1;
AdTrue4(7,5) = 1;  AdTrue4(5,7) = 1;
AdTrue4(7,3) = 1;  AdTrue4(3,7) = 1;
AdTrue4(7,12) = 1;  AdTrue4(12,7) = 1;
AdTrue4(5,1) = 1;  AdTrue4(1,5) = 1;
AdTrue4(1,17) = 1;  AdTrue4(17,1) = 1;
AdTrue4(17,15) = 1;  AdTrue4(15,17) = 1;
AdTrue4(15,2) = 1;  AdTrue4(2,15) = 1;
AdTrue4(13,2) = 1;  AdTrue4(2,13) = 1;
AdTrue4(15,14) = 1;  AdTrue4(14,15) = 1;
AdTrue4(14,2) = 1;  AdTrue4(2,14) = 1;
AdTrue4(17,4) = 1;  AdTrue4(4,17) = 1;
AdTrue4(2,19) = 1;  AdTrue4(19,2) = 1;
AdTrue4(4,11) = 1;  AdTrue4(11,4) = 1;
AdTrue4(4,16) = 1;  AdTrue4(16,4) = 1;
AdTrue4(19,3) = 1;  AdTrue4(3,19) = 1;
AdTrue4(3,12) = 1;  AdTrue4(12,3) = 1;
AdTrue4(3,18) = 1;  AdTrue4(3,18) = 1;
AdTrue4(18,20) = 1;  AdTrue4(20,18) = 1;

%hyperpar 2
lambda1G=0.25;
lambda2G=60;
display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

Ad1 = ZG(:,:,1)~=0;
Ad2 = ZG(:,:,60)~=0;
Ad3 = ZG(:,:,100)~=0;
Ad4 = ZG(:,:,150)~=0;

TP = TP + sum(sum(Ad1 .* AdTrue1))/2 + sum(sum(Ad2 .* AdTrue2))/2 + sum(sum(Ad3 .* AdTrue3))/2 + sum(sum(Ad4 .* AdTrue4))/2;
 FN = FN + sum(sum((1-Ad1) .* AdTrue1))/2 + sum(sum((1-Ad2) .* AdTrue2))/2 + sum(sum((1-Ad3) .* AdTrue3))/2 + sum(sum((1-Ad4) .* AdTrue4))/2;
 TN = TN + (sum(sum((1-Ad1) .* (1-AdTrue1)))-20 + sum(sum((1-Ad2) .* (1-AdTrue2)))-20 + sum(sum((1-Ad3) .* (1-AdTrue3)))-20 + sum(sum((1-Ad4) .* (1-AdTrue4)))-20 )/2;
 FP = FP + sum(sum((Ad1) .* (1-AdTrue1)))/2 + sum(sum((Ad2) .* (1-AdTrue2)))/2 + sum(sum((Ad3) .* (1-AdTrue3)))/2 + sum(sum((Ad4) .* (1-AdTrue4)))/2 ;

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_scenario4_0.2.png','Resolution',300)

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,4,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,4,2);
plotGraph(squeeze(abs(ZG(:,:,(60)))),9,0);
subplot(1,4,3);
plotGraph(squeeze(abs(ZG(:,:,(100)))),9,0);
subplot(1,4,4);
plotGraph(squeeze(abs(ZG(:,:,(150)))),9,0);
exportgraphics(gcf,'G_scenario4_0.2.png','Resolution',300)

%hyperpar 3
lambda1G=0.1;
lambda2G=60;
display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_scenario4_0.1.png','Resolution',300)

%% Graphical estimation

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,4,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,4,2);
plotGraph(squeeze(abs(ZG(:,:,(60)))),9,0);
subplot(1,4,3);
plotGraph(squeeze(abs(ZG(:,:,(63)))),9,0);
subplot(1,4,4);
plotGraph(squeeze(abs(ZG(:,:,(100)))),9,0);
exportgraphics(gcf,'G_scenario4_0.1.png','Resolution',300)

%% Read Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = csvread("C:\Users\beatr\Documents\weekly.csv")

lambda1G=0.25;
lambda2G=20;
%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.25_20.png','Resolution',300)

%% Graphical estimation

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,4,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,4,2);
plotGraph(squeeze(abs(ZG(:,:,(61)))),9,0);
subplot(1,4,3);
plotGraph(squeeze(abs(ZG(:,:,(80)))),9,0);
subplot(1,4,4);
plotGraph(squeeze(abs(ZG(:,:,(98)))),9,0);
exportgraphics(gcf,'G_0.25_20.png','Resolution',300)

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,80))
clustCoeff_modified(ZG(:,:,98))

%n.2
lambda1G=0.35;
lambda2G=20;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.35_20.png','Resolution',300)
%% Graphical estimation

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,5,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,5,2);
plotGraph(squeeze(abs(ZG(:,:,(61)))),9,0);
subplot(1,5,3);
plotGraph(squeeze(abs(ZG(:,:,(80)))),9,0);
subplot(1,5,4);
plotGraph(squeeze(abs(ZG(:,:,(98)))),9,0);
subplot(1,5,5);
plotGraph(squeeze(abs(ZG(:,:,(99)))),9,0);
exportgraphics(gcf,'G_0.35_20.png','Resolution',300)

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,80))
clustCoeff_modified(ZG(:,:,98))
clustCoeff_modified(ZG(:,:,99))

%n.3
lambda1G=0.5;
lambda2G=20;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.5_20.png','Resolution',300)

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,5,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,5,2);
plotGraph(squeeze(abs(ZG(:,:,(61)))),9,0);
subplot(1,5,3);
plotGraph(squeeze(abs(ZG(:,:,(68)))),9,0);
subplot(1,5,4);
plotGraph(squeeze(abs(ZG(:,:,(77)))),9,0);
subplot(1,5,5);
plotGraph(squeeze(abs(ZG(:,:,(80)))),9,0);
exportgraphics(gcf,'G_0.5_20.png','Resolution',300)

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,68))
clustCoeff_modified(ZG(:,:,77))
clustCoeff_modified(ZG(:,:,80))

%n.4
lambda1G=0.25;
lambda2G=10;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.25_10.png','Resolution',300) 

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,57))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,80))
clustCoeff_modified(ZG(:,:,98))
clustCoeff_modified(ZG(:,:,99))
clustCoeff_modified(ZG(:,:,116))

%n.5
lambda1G=0.35;
lambda2G=10;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.35_10.png','Resolution',300) 

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,5,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,5,2);
plotGraph(squeeze(abs(ZG(:,:,(61)))),9,0);
subplot(1,5,3);
plotGraph(squeeze(abs(ZG(:,:,(80)))),9,0);
subplot(1,5,4);
plotGraph(squeeze(abs(ZG(:,:,(98)))),9,0);
subplot(1,5,5);
plotGraph(squeeze(abs(ZG(:,:,(116)))),9,0);
exportgraphics(gcf,'G_0.35_10.png','Resolution',300)

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,80))
clustCoeff_modified(ZG(:,:,98))
clustCoeff_modified(ZG(:,:,116))

%n.6
lambda1G=0.5;
lambda2G=10;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.5_10.png','Resolution',300) 

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,57))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,69))
clustCoeff_modified(ZG(:,:,70))
clustCoeff_modified(ZG(:,:,78))
clustCoeff_modified(ZG(:,:,80))

%n.7
lambda1G=0.25;
lambda2G=60;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.25_60.png','Resolution',300) 

clustCoeff_modified(ZG(:,:,1))


%n.8
lambda1G=0.35;
lambda2G=60;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.35_60.png','Resolution',300) 

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,3,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,3,2);
plotGraph(squeeze(abs(ZG(:,:,(61)))),9,0);
subplot(1,3,3);
plotGraph(squeeze(abs(ZG(:,:,(80)))),9,0)
exportgraphics(gcf,'G_0.35_60.png','Resolution',300)

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,80))

%n.9
lambda1G=0.5;
lambda2G=60;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.5_60.png','Resolution',300) 

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,3,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,3,2);
plotGraph(squeeze(abs(ZG(:,:,(61)))),9,0);
subplot(1,3,3);
plotGraph(squeeze(abs(ZG(:,:,(80)))),9,0)
exportgraphics(gcf,'G_0.5_60.png','Resolution',300)

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,80))

%n.10
lambda1G=0.25;
lambda2G=55;

%% Read Data
%y = csvread("C:\Users\beatr\Documents\weekly.csv")

%% RUN GFGL
% As discussed in the paper, our implementation of GFGL is highly dependent
% on the number of changepoints (K) estimated. Increasing lambda2G
% decreases the number of estimated changepoints.

display('Running GFGL');
tic
[ ThetaG,ZG,cpG,SG,~,~ ] = GFGL( y,lambda1G,lambda2G, 1);
disp(cpG)
t=toc;
display(['GFGL took t=',int2str(t),'seconds to find K=',num2str(length(cpG))]);

%% Output Results

%% Changepoints and Edge estimates

figure(1)
set(gcf, 'Color', [1,1,1]);
cpI = []
plotCP(cpG+1,cpI,y,0);
title('change points detected with GFGL');
exportgraphics(gcf,'cp_0.25_55.png','Resolution',300) 

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(1,3,1);
plotGraph(squeeze(abs(ZG(:,:,1))),9,0);
title('GFGL estimate');
subplot(1,3,2);
plotGraph(squeeze(abs(ZG(:,:,(61)))),9,0);
subplot(1,3,3);
plotGraph(squeeze(abs(ZG(:,:,(80)))),9,0)
exportgraphics(gcf,'G_0.25_55.png','Resolution',300)

clustCoeff_modified(ZG(:,:,1))
clustCoeff_modified(ZG(:,:,61))
clustCoeff_modified(ZG(:,:,80))