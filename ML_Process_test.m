% Organizing some ML stuff
% Author: Jonathan Weinrib
% Originated 01/21/2019

%% This file has some work starting with the ML process once I already have
% The dataset fully formed
clc, clear all, close all

%%
load('signal_list.mat')
sig_1_orig = load('est_enf_B_P1.mat');
sig_2_orig = load('est_enf_B_P2.mat');
sig_3_orig = load('est_enf_D_P1.mat');
sig_4_orig = load('est_enf_D_P2.mat');
%%


x1 = processed_ENF_cellArray{1}{4};
x1labels = [ones(length(x1),1) ones(length(x1),1)];

x2 = processed_ENF_cellArray{2}{4};
x2labels = [ones(length(x2),1) ones(length(x2),1)];

x3 = processed_ENF_cellArray{3}{4};
x3labels = [ones(length(x3),1) ones(length(x3),1)];


x4 = processed_ENF_cellArray{4}{4};
x4labels = [ones(length(x4),1) ones(length(x4),1)];


x5 = processed_noise_cellArray{1}{4};
x5labels = [zeros(length(x5),1) ones(length(x5),1)];

x6 = processed_noise_cellArray{2}{4};
x6labels = [zeros(length(x6),1) ones(length(x6),1)];
x7 = processed_noise_cellArray{3}{4};
x7labels = [zeros(length(x7),1) ones(length(x7),1)];
x8 = processed_noise_cellArray{4}{4};
x8labels = [zeros(length(x8),1) ones(length(x8),1)];







%Labels:: % first column is ENF, second column is noise
YtraiN = {categorical(x1labels(:,1)') categorical(x3labels(:,1)') categorical(x5labels(:,1)') categorical(x7labels(:,1)')};
YtesT = {categorical(x2labels(:,1)') categorical(x4labels(:,1)') categorical(x6labels(:,1)') categorical(x8labels(:,1)')};



xtrain = {x1 x3 x5 x7};
xtest = {x2 x4 x6 x8};


%% make training and testing data
% Time to make xtrain and ytrain




%% Logistic Regression Learner:
%y is a 4x2 matrix
x1labels = [1 1];
x2labels = [1 1];
x3labels = [1 1];
x4labels = [1 1];
x5labels = [0 1];
x6labels = [0 1];
x7labels = [0 1];
x8labels = [0 1];

ytrain = [x1labels;x3labels;x5labels;x7labels];
ytest = [x2labels;x4labels;x6labels;x8labels];
%%
% take global features from each of the estimated ENF's
% x should be a numrecordingsxnumfeatures matrix

xtrain = {x1 x3 x5 x7};
xtest = {x2 x4 x6 x8};

Xtrain = [];
Xtest = [];
for k = 1:length(xtrain)
    Xtrain(k,1) = mean(xtrain{k});
    Xtrain(k,2) = var(xtrain{k});
    
    Xtest(k,1) = mean(xtest{k});
    Xtest(k,2) = var(xtest{k});
    
end

%%
% Now that I have a basic feature matrix, let's plug into classifier:
% Logistic Regression
Log_Reg_Coeff = mnrfit(Xtrain,ytrain);
logReg_Preds = mnrval(Log_Reg_Coeff,Xtest);

%% Discrim Analysis:
%lda_0 = fitcdiscr(Xtrain,ytrain);


%%
%ar_model_0 = ar(

%% NNets
featureDim = 1;
numHiddenUnits = 100;
numClasses = 2;

layers = [ ...
    sequenceInputLayer(featureDim)
    lstmLayer(numHiddenUnits,'OutputMode','sequence')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];


options = trainingOptions('adam', ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',20, ...
    'Verbose',0, ...
    'Plots','training-progress');

%ytrain_nnet = ytrain(:,1);
%ytest_nnet = ytest(:,1);

net_0 = trainNetwork(xtrain,YtraiN,layers,options);

%% Plotting and viewing

% x = xtrain{1};
% %classes = 
% X = XTrain{1}(1,:);
% classes = categories(YTrain{1});
% 
% figure
% for j = 1:numel(classes)
%     label = classes(j);
%     idx = find(YTrain{1} == label);
%     hold on
%     plot(idx,X(idx))
% end
% hold off
% 
% xlabel("Time Step")
% ylabel("Acceleration")
% title("Training Sequence 1, Feature 1")
% legend(classes,'Location','northwest')

