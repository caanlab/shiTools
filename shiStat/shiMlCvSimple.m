function [Error,Info] = shiMlCvSimple(Data,funMdl,funHat,funErr,KFold,McRep,StratVar,verbose)

% performes machine learning with simple cross-validation (no hyperparameter optimization)
%
% Input:
%      Data     : input Y and X data together, letting funMdl and funHat decide which columns are Y and which are X (so as to handle unsupervised learning with X only)
%      funMdl   : handle for training function of form       Mdl=func(Data)     , run once in each CV partition, where Data [ndims=2, e.g. nobs-by-nvar matrix] contains training data, and Mdl contains the trained predictive model parameters
%      funHat   : handle for testing function of form        Hat=func(Data,Mdl) , run once in each CV partition, where Data [ndims=2, e.g. nobs-by-nvar matrix] contains testing data, Mdl is obtained from funMdl, and Hat contains any information about the model and its performance [which will be packed in a cell and passed to funErr]
%      funErr   : handle for error summary function of form  Err=func(HatCell)  , to summarize predictive performance across CV partitions, where each element of HatCell [ncv-by-1 cell] contains Hat from funHat on one CV partition, and Err is a scalar aggregate error value such as MSE or 1-ROC_AUC (default = @(x)mean([x{:}]))
%      KFold    : k-fold cross-validation for prediction (default = 10)
%      McRep    : Monte-Carlo repetitions of k-fold partitioning for prediction (default = 1, i.e. k-fold only once, no further Monte-Carlo)
%      StratVar : stratification variable (see cvpartition.m; default = [])
%
% Output:
%      Error    : scalar, summarized prediction error across CV partitions returned by funErr
%      Info     : structure containing the following information from each CV partition
%                 .McIdx         - iteration index of Monte-Carlo repetitions, size = [KFold*McRep, 1]
%                 .CvIdx         - iteration index of k-fold CV partitions within each Monte-Carlo repitition, size = [KFold*McRep, 1]
%                 .Hat           - cell array of model/prediction information for each CV partition returned by funHat, size = [KFold*McRep, 1]
%
% Example:
%     Data = [randn(50,1),randn(50,6)]; % 1 Y variable, 6 X variables, 50 observations
%     funMdl = @(DATA)fitglm(DATA(:,2:6),DATA(:,1)); % model training function (e.g. linear regression) with data as input, returning a model variable that's readable by funHat below
%     funHat = @(DATA,MDL)[DATA(:,1),predict(MDL,DATA(:,2:6))]; % model testing function with data and trained model (from funMdl) as inputs, returning any information that needs to be returned and is requred for error calculation by funErr (e.g. both y and y_hat)
%     funErr = @(HATCELL)mean(diff(cat(1,HATCELL{:}),[],2).^2); % prediction error calculation function whose input is a cell array (each cell contains an output from funHat), returning a scalar prediction error (e.g. MSE)
%     KFold = 50; % e.g. leave-one-out CV
%     McRep =  1; % e.g. no Monte-Carlo for CV
%     StratVar = randn(50,1)>0; % optional, e.g. stratify CV partitions by a binary (possibly gender) variable
%     [Error,Info] = shiMlNested(Data,funMdl,funHat,funErr,KFold,McRep,StratVar)
%
% Zhenhao Shi 2020-10-03

%% initializing input

if ~exist('funErr','var') || isempty(funErr)
    funErr = @(x)mean([x{:}]); % default: average all
end

if ~exist('KFold','var') || isempty(KFold)
    KFold = 10; % default: 10-fold CV for prediction
end

if ~exist('McRep','var') || isempty(McRep)
    McRep = 1; % default: 1 Monte-Carlo iteration of CV for prediction
end

if ~exist('StratVar','var') || isempty(StratVar)
    Group = size(Data,1); % default: left empty, do not stratify during CVs
else
    Group = StratVar;
end

if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

%% defining CV partitions

if isinf(KFold)
    KFold = size(Data,1);
    McRep = 1;
end

[isTrain,idxMc,idxCv] = shiMlCvPartition(Group,KFold,McRep);
Info.McIdx = idxMc';
Info.CvIdx = idxCv';

%% main job

Info.Hat = cell(KFold*McRep, 1);

if verbose, fprintf('Running simple cross-validation...   %%'); end
for iter = 1:(KFold*McRep) % loop through all prediction CV partitions
    if verbose, fprintf('\b\b\b\b%3d%%',round(iter/(KFold*McRep)*100)); end
    
    %% creating training and testing datasets

    DataTrain = Data(isTrain(:,iter),:);
    DataTest = Data(~isTrain(:,iter),:);

    %% prediction
    
    Mdl = funMdl(DataTrain);
    Hat = funHat(DataTest,Mdl);

    Info.Hat{iter,1} = Hat;
    
end

%% final prediction error aggregated across CV partitions

Error = funErr(Info.Hat);
if verbose, display(Error); end


function [isTrain,idxMc,idxCv] = shiMlCvPartition(Group,KFold,McRep)
%%
N = shiIf(isscalar(Group),Group,length(Group));
isTrain = false(N,KFold*McRep);
idxMc = nan(1,KFold*McRep);
idxCv = nan(1,KFold*McRep);
for mc = 1:McRep
    rng shuffle;
    CV = cvpartition(Group,'KFold',KFold);
    for cv = 1:KFold
        idxMc((mc-1)*KFold+cv) = mc;
        idxCv((mc-1)*KFold+cv) = cv;
        isTrain(:,(mc-1)*KFold+cv) = training(CV,cv);
    end
end