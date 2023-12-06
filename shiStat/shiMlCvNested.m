function [Error,Info,InfoAll] = shiMlCvNested(Data,funMdl,funHat,funErr,ParamSet,KFold_Outer,McRep_Outer,KFold_Inner,McRep_Inner,StopOptimAtRise,StratVar)

% (beta) performes machine learning with nested cross-validation for hyperparameter optimization
%
% Input:
%      Data            : input Y and X data together, letting funMdl and funHat to decide which columns are Y and which are X (so as to handle unsupervised learning with X only)
%      funMdl          : handle for training function of form       Mdl=func(Data,Param) , run once in each CV partition, where Data [ndims=2, e.g. nobs-by-nvar matrix] contains training data, Param contains parameters to be optimized, and Mdl contains the trained predictive model parameters
%      funHat          : handle for testing function of form        Hat=func(Data,Mdl)   , run once in each CV partition, where Data [ndims=2, e.g. nobs-by-nvar matrix] contains testing data, Mdl is obtained from funMdl, and Hat contains any information about the model and its performance [which will be packed in a cell and passed to funErr]
%      funErr          : handle for error summary function of form  Err=func(HatCell)    , to summarize predictive performance across CV partitions, where each element of HatCell [ncv-by-1 cell] contains Hat from funHat on one CV partition, and Err is a scalar aggregate error value such as MSE or 1-ROC_AUC (default = @(x)mean([x{:}]))
%      ParamSet        : cell matrix for parameters used by funMdl, each row specifying a unique combination of parameters, typically resulting in increasing model complexity, letting nested CV to determine the optimal set (i.e. hyperparameter optimization)
%      KFold_Outer     : k-fold cross-validation for prediction (outer loop) (default = 10)
%      McRep_Outer     : Monte-Carlo repetitions of k-fold partitioning for prediction (outer loop) (default = 1, i.e. k-fold only once, no further Monte-Carlo)
%      KFold_Inner     : k-fold cross-validation for hyperparameter optimization (inner loop) (default = 10)
%      McRep_Inner     : Monte-Carlo repetitions of k-fold partitioning for hyperparameter optimization (inner loop) (default = 1, i.e. k-fold only once, no further Monte-Carlo)
%      StopOptimAtRise : whether to stop evaluating further parameter sets if prediction error (from funErr) starts to increase (default = false, i.e. greedy testing of all parameters)
%      StratVar        : stratification variable (see cvpartition.m; default = [])
%
% Output:
%      Error           : scalar, summarized prediction error across outer CV partitions returned by funErr
%      Info            : structure containing the following information from each outer CV partition
%                         .McIdx         - iteration index of outer Monte-Carlo repetitions, size = [KFold_Outer*McRep_Outer, 1]
%                         .CvIdx         - iteration index of outer k-fold CV partitions within each Monte-Carlo repitition, size = [KFold_Outer*McRep_Outer, 1]
%                         .Error_InnerCv - prediction error for each outer CV partition and each parameter set, size = [KFold_Outer*McRep_Outer, size(ParamSet,1)]
%                         .OptParam      - index of best parameter set used in each outer CV partition, size = [KFold_Outer*McRep_Outer, 1]
%                         .Hat           - cell array of model/prediction information for each outer CV partition returned by funHat, size = [KFold_Outer*McRep_Outer, 1]
%      InfoAll         : structure containing the following information from each inner CV partition
%                         .OuterMcIdx    - iteration index of outer Monte-Carlo repetitions, size = [KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1]
%                         .OuterCvIdx    - iteration index of outer k-fold CV partitions, size = [KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1]
%                         .InnerMcIdx    - iteration index of inner Monte-Carlo repetitions, size = [KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1]
%                         .InnerCvIdx    - iteration index of inner k-fold CV partitions, size = [KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1]
%                         .Hat           - cell array of model/prediction information for each inner CV partition returned by funHat, size = [KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1]
%
% Example:
%     Data = [randn(50,1),randn(50,6)]; % 1 Y variable, 6 X variables, 50 observations
%     funMdl = @(DATA,PARAM)fitglm(DATA(:,2:(1+PARAM)),DATA(:,1)); % model training function (e.g. linear regression) with data and parameter (e.g. number of predictors) as inputs, returning a model variable that's readable by funHat below
%     funHat = @(DATA,MDL)[DATA(:,1),predict(MDL,DATA(:,2:(1+MDL.NumPredictors)))]; % model testing function with data and trained model (from funMdl) as inputs, returning any information that needs to be returned and is requred for error calculation by funErr (e.g. both y and y_hat)
%     funErr = @(HATCELL)mean(diff(cat(1,HATCELL{:}),[],2).^2); % prediction error calculation function whose input is a cell array (each cell contains an output from funHat), returning a scalar prediction error (e.g. MSE)
%     ParamSet = {1;2;3;4;5;6}; % possible parameters for funMdl (e.g. number of predictors), which will be tuned
%     KFold_Outer = 50; % e.g. leave-one-out for outer CV
%     McRep_Outer =  1; % e.g. no Monte-Carlo for outer CV
%     KFold_Inner = 10; % e.g. 10-fold for inner CV
%     McRep_Inner = 50; % e.g. 50 Monte-Carlo repetitions of 10-fold for inner CV, thus 500 partitions in total
%     StopOptimAtRise = true; % e.g. stop testing more complex linear regression models when prediction error starts increasing
%     StratVar = randn(50,1)>0; % optional, e.g. stratify CV partitions by a binary (possibly gender) variable
%     [Error,Info,InfoAll] = shiMlNested(Data,funMdl,funHat,funErr,ParamSet,KFold_Outer,McRep_Outer,KFold_Inner,McRep_Inner,StopOptimAtRise,StratVar)
%
% Zhenhao Shi 2020-10-03

%% initializing input

if ~exist('funErr','var') || isempty(funErr)
    funErr = @(x)mean([x{:}]); % default: average all
end

if ~exist('KFold_Outer','var') || isempty(KFold_Outer)
    KFold_Outer = 10; % default: 10-fold CV for prediction
end

if ~exist('McRep_Outer','var') || isempty(McRep_Outer)
    McRep_Outer = 1; % default: 1 Monte-Carlo iteration of CV for prediction
end

if ~exist('KFold_Inner','var') || isempty(KFold_Inner)
    KFold_Inner = 10; % default: 10-fold CV for hyperparameter optimization
end

if ~exist('McRep_Inner','var') || isempty(McRep_Inner)
    McRep_Inner = 1; % default: 1 Monte-Carlo iteration of CV for hyperparameter optimization
end

if ~exist('StopOptimAtRise','var') || isempty(StopOptimAtRise)
    StopOptimAtRise = false; % default: traverse the entire hyperparameter set without stopping at rising error
end

if ~exist('StratVar','var') || isempty(StratVar)
    Group = size(Data,1); % default: left empty, do not stratify during CVs
else
    Group = StratVar;
end

%% checking if output #3 would be too large

if nargin==3
    xxtrain = 1:round(size(Data,1)*(KFold_Outer-1)/KFold_Outer);
    xxtest = setdiff(1:size(Data,1),xxtrain);
    xxmdl = funMdl(Data(xxtrain,:),ParamSet(end,:)); % try last parameter set
    xxhat = {funHat(Data(xxtest,:),xxmdl)}; %#ok<NASGU>
    xx = getfield(whos('xxhat'),'bytes');
    xx = xx*KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner/1000000000;
    if xx>5 || xx<15
        warning('output #3 from shiMlNested will be over %.2f GB',xx);
    elseif xx>=15
        shiDisp(sprintf('output #3 from shiMlNested will be over %.2f GB',xx),true);
    end
end

%% defining outer CV partitions

[isTrain_Outer,idxMc_Outer,idxCv_Outer] = shiMlCvPartition(Group,KFold_Outer,McRep_Outer);

%% main job

if nargout >= 3
    InfoAll.OuterMcIdx = nan(KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1);
    InfoAll.OuterCvIdx = nan(KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1);
    InfoAll.InnerMcIdx = nan(KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1);
    InfoAll.InnerCvIdx = nan(KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, 1);
    InfoAll.Hat = cell(KFold_Outer*McRep_Outer*KFold_Inner*McRep_Inner, size(ParamSet,1));
end

Info.McIdx = nan(KFold_Outer*McRep_Outer, 1);
Info.CvIdx = nan(KFold_Outer*McRep_Outer, 1);
Info.OptParam = nan(KFold_Outer*McRep_Outer, 1);
Info.Error_InnerCv = nan(KFold_Outer*McRep_Outer, size(ParamSet,1));
Info.Hat = cell(KFold_Outer*McRep_Outer, 1);

fprintf('Running nested cross-validation...   %%');
for iterOut = 1:(KFold_Outer*McRep_Outer) % loop through all prediction CV partitions
    fprintf('\b\b\b\b%3d%%',round(iterOut/(KFold_Outer*McRep_Outer)*100));
    
    %% creating training and testing datasets

    DataTrain = Data(isTrain_Outer(:,iterOut),:);
    DataTest = Data(~isTrain_Outer(:,iterOut),:);

    if isscalar(Group)
        GroupTrain = sum(isTrain_Outer(:,iterOut));
    else
        GroupTrain = Group(isTrain_Outer(:,iterOut));
    end

    %% defining inner CV partitions

    [isTrain_Inner,idxMc_Inner,idxCv_Inner] = shiMlCvPartition(GroupTrain,KFold_Inner,McRep_Inner);
 
    if nargout >= 3
        InfoAll.OuterMcIdx(( (iterOut-1)*(KFold_Inner*McRep_Inner)+1 ) : ( iterOut*(KFold_Inner*McRep_Inner) )) = idxMc_Outer(iterOut);
        InfoAll.OuterCvIdx(( (iterOut-1)*(KFold_Inner*McRep_Inner)+1 ) : ( iterOut*(KFold_Inner*McRep_Inner) )) = idxCv_Outer(iterOut);
        InfoAll.InnerMcIdx(( (iterOut-1)*(KFold_Inner*McRep_Inner)+1 ) : ( iterOut*(KFold_Inner*McRep_Inner) )) = idxMc_Inner;
        InfoAll.InnerCvIdx(( (iterOut-1)*(KFold_Inner*McRep_Inner)+1 ) : ( iterOut*(KFold_Inner*McRep_Inner) )) = idxCv_Inner;
    end
    
    Info.McIdx(iterOut) = idxMc_Outer(iterOut);
    Info.CvIdx(iterOut) = idxCv_Outer(iterOut);
    
    %% searching through parameter sets
    
    for p = 1:size(ParamSet,1) % loop through all possible parameters on the mcv-th outer CV partition

        xHat = [];
        for iterIn = 1:(KFold_Inner*McRep_Inner) % loop through all optimization CV partitions to test the p-th parameter set
            xMdl = funMdl(DataTrain(isTrain_Inner(:,iterIn),:),ParamSet(p,:));
            xHat{iterIn,1} = funHat(DataTrain(~isTrain_Inner(:,iterIn),:),xMdl); %#ok<*AGROW> % store hat information in a cell
            if nargout >= 3
                InfoAll.Hat((iterOut-1)*(KFold_Inner*McRep_Inner)+iterIn,p) = xHat(iterIn);
            end
        end % done getting hats for p-th parameter set, error information for each inner CV partition stored in xErr

        Info.Error_InnerCv(iterOut,p) = funErr(xHat);

        if StopOptimAtRise && p>1 && Info.Error_InnerCv(iterOut,p)>Info.Error_InnerCv(iterOut,p-1)
            break;
        end

    end % done getting aggregated errors for all parameter sets, error value stored in InfoOut.Error_InnerCv(iterOut,:)
    
    %% prediction using the best parameter set
    
    [~, Info.OptParam(iterOut)] = nanmin(Info.Error_InnerCv(iterOut,:));
    Mdl = funMdl(DataTrain,ParamSet(Info.OptParam(iterOut),:));
    Hat = funHat(DataTest,Mdl);

    Info.Hat{iterOut,1} = Hat;
    
end

%% final prediction error aggregated across outer CV partitions

Error = funErr(Info.Hat);
fprintf('.  Error = %g\n',Error);


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