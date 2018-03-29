
%% 
% this cross validation script is only used for 'forward' and 'backward'
% states.

%% clean working space
% clear all;
% close all;
% clc;

%% load data
% before this part, make sure all files exist
% if 'dataX'/'dataY' is missing, run 'genTraining.m' first
% load alldata;
% load dataX;
% load dataY;

%% initialize parameters
% expNum = size(alldata, 2);
% frameLen = 100;
% posLabel = 'fwd';
% startOffset = 0;
% endOffset = 0;
% channels = size(alldata(1).dFF_smooth, 1);
% CVFolds = 5;
% CVResults = [];

%% set split
% setIndxs = CVSplit(expNum, CVFolds);

%% start loop
for i = 1:CVFolds
%     fprintf(['cross validation ' num2str(i) ' is running ...\n']);
%     trainingIndxs = [];
%     testingIndxs = setIndxs{i};
%     for j = 1:CVFolds
%         if(j==i)
%             continue;
%         end
%         trainingIndxs = [trainingIndxs setIndxs{j}]; 
%     end
    [trainingX, trainingY, testingY, frameId, expId] = genTrainingData(alldata, dataX, dataY, trainingIndxs, testingIndxs, posLabel, startOffset, endOffset);
    fprintf(['experiment set ' num2str(sort(trainingIndxs)) ' are used for training.\n']);
    fprintf(['experiment set ' num2str(sort(testingIndxs)) ' are used for testing.\n']);
    trainingMean = mean(trainingX, 1);
    trainingStd = std(trainingX, 0, 1);
    svmModels = svmTrain((trainingX-trainingMean)./trainingStd, trainingY);
    prediction = svmPredict(alldata, testingIndxs, svmModels, posLabel, frameLen, trainingMean, trainingStd);
    [CVStats, mismatchFrame, mismatchExp, mismatchIndx] = computeTestResults(prediction, testingY, frameId, expId);
    pause;
end

%% functions added

% ------------------------------------------------------------------------------- %
% function CVSplit
%       randomly split the experiments to k folds to perform k-fold
%       validation
%       this function does not split the training set or testing set, it
%       only gives the indexs for splitting
% ------------------------------------------------------------------------------- % 
function setIndxs = CVSplit(expNum, CVFolds)
    setNum = floor(expNum/CVFolds);
    indxs = linspace(1, expNum, expNum);
    setIndxs = cell(CVFolds, 1);
    for i = 1:CVFolds-1
        setIndxs{i} = randsample(indxs, setNum);
        indxs = setdiff(indxs, setIndxs{i});
    end
    setIndxs{CVFolds} = indxs(randperm(length(indxs)));
end

% ------------------------------------------------------------------------------- %
% function genTrainingData
%       generate training data & label, testing data & label based on
%       splitted indexs
% ------------------------------------------------------------------------------- % 
function [trainingX, trainingY, testingY, frameId, expId] = genTrainingData(alldata, dataX, dataY, trainingIndxs, testingIndxs, posLabel, startOffset, endOffset)
    trainingX = [];
    trainingY = [];
    testingY = [];
    
    for i = 1:length(trainingIndxs)
        trainingX = [trainingX; dataX{trainingIndxs(i)}];
        trainingY = [trainingY; dataY{trainingIndxs(i)}];
    end
    [testingY, frameId, expId] = genLabels(alldata, testingIndxs, posLabel, startOffset, endOffset);
end

% ------------------------------------------------------------------------------- %
% function genLabels
% ------------------------------------------------------------------------------- % 
function [labels, frameId, expId] = genLabels(alldata, testingIndxs, posLabel, startOffset, endOffset)
    assert(startOffset>=0 && endOffset>=0);
    labels = [];
    frameId = [];
    expId = [];
    switch(posLabel)
        case 'fwd'
            posLabel = 1;
        case 'bkd'
            posLabel = 3;
        otherwise
            posLabel = -1;
    end
    assert(posLabel~=-1);
    for i = 1:length(testingIndxs)
        curLabels = zeros(length(alldata(testingIndxs(i)).dFF_smooth), 1);
        tempFrameId = linspace(1, length(curLabels), length(curLabels));
        tempExpId = ones(size(curLabels))*testingIndxs(i);
        annos = alldata(testingIndxs(i)).anno;
        for j = 1:length(annos)
            startIndx = annos(j, 1)+1;
            endIndx = annos(j, 2)+1;
            curLabel = annos(j, 3);
            if(curLabel==posLabel)
                label = 1;
            elseif(curLabel==2 && posLabel==1)
                label = NaN;
            elseif(curLabel==4 && posLabel==3)
                label = NaN;
            else
                continue;
            end
            curLabels(startIndx:endIndx) = label;
            curLabels(startIndx:startIndx+startOffset) = NaN;
            curLabels(endIndx-endOffset:endIndx) = NaN;
        end
        labels = [labels; curLabels];
        frameId = [frameId; tempFrameId'];
        expId = [expId; tempExpId];
    end
end

% ------------------------------------------------------------------------------- % 
% function trainSVM
%       call fitcsvm to train classifier, 
%       methods for processing like 'balanced' are not included for now
% ------------------------------------------------------------------------------- % 
function SVMModel = svmTrain(trainingX, trainingY)
    countY = length(trainingY);
    count0 = sum(trainingY==0);
    count1 = sum(trainingY==1);
    assert(count0+count1 == countY);
    
    weightVec = zeros(size(trainingY));
    weightVec(trainingY==0) = countY / count0 / 2;
    weightVec(trainingY==1) = countY / count1 / 2;
               
    SVMModel = fitcsvm(trainingX, trainingY, 'Weights', weightVec, 'cost', [0 1; 1 0]);
end

% ------------------------------------------------------------------------------- % 
% function computeTestResults
%       evaluate performance of classifier
% ------------------------------------------------------------------------------- % 
function [CVStats, mismatchFrame, mismatchExp, mismatchIndx] = computeTestResults(prediction, testingY, frameId, expId)
    mismatchFrame = [];
    mismatchExp = [];
    mismatchIndx = [];
    pos_pos = 0;
    pos_neg = 0;
    pos_ref = 0;
    neg_pos = 0;
    neg_neg = 0;
    neg_ref = 0;
    nan_pos = 0;
    nan_neg = 0;
    nan_ref = 0;
    for i = 1:length(testingY)
        if(isnan(testingY(i)))
            if(isnan(prediction(i)))
                nan_ref = nan_ref + 1;
            elseif(prediction(i)==1)
                nan_pos = nan_pos + 1;
                mismatchIndx = [mismatchIndx i];
                mismatchFrame = [mismatchFrame frameId(i)];
                mismatchExp = [mismatchExp expId(i)];
            else
                nan_neg = nan_neg + 1;
                mismatchIndx = [mismatchIndx i];
                mismatchFrame = [mismatchFrame frameId(i)];
                mismatchExp = [mismatchExp expId(i)];
            end
        elseif(testingY(i)==1)
            if(isnan(prediction(i)))
                pos_ref = pos_ref + 1;
                mismatchIndx = [mismatchIndx i];
                mismatchFrame = [mismatchFrame frameId(i)];
                mismatchExp = [mismatchExp expId(i)];
            elseif(prediction(i)==1)
                pos_pos = pos_pos + 1;
            else
                pos_neg = pos_neg + 1;
                mismatchIndx = [mismatchIndx i];
                mismatchFrame = [mismatchFrame frameId(i)];
                mismatchExp = [mismatchExp expId(i)];
            end
        else
            if(isnan(prediction(i)))
                neg_ref = neg_ref + 1;
                mismatchIndx = [mismatchIndx i];
                mismatchFrame = [mismatchFrame frameId(i)];
                mismatchExp = [mismatchExp expId(i)];
            elseif(prediction(i)==1)
                neg_pos = neg_pos + 1;
                mismatchIndx = [mismatchIndx i];
                mismatchFrame = [mismatchFrame frameId(i)];
                mismatchExp = [mismatchExp expId(i)];
            else
                neg_neg = neg_neg + 1;
            end
        end
    end
    CVStats = [pos_pos, pos_neg, pos_ref, neg_pos, neg_neg, neg_ref, nan_pos, nan_neg, nan_ref];
end

% ------------------------------------------------------------------------------- % 
% function svmPredict
%       make prediction using svm classifier
% ------------------------------------------------------------------------------- % 
function prediction = svmPredict(alldata, testingIndxs, svmModels, posLabel, frameLen, trainingMean, trainingStd)
    assert(prod(posLabel=='fwd') || prod(posLabel=='bkd'));
    prediction = [];
    for i = 1:length(testingIndxs)
        curPrediction = nan(length(alldata(testingIndxs(i)).dFF_smooth), 1);
        data = alldata(testingIndxs(i)).dFF_smooth;
        for j = 1:length(data)
            if(~quick_pred_active(data, j))
                curPrediction(j) = 0;
                continue;
            end
            if(j<frameLen)
                continue;
            end
            fval_0 = compute_special_feature(data, j);
            fval_1 = compute_special_feature(data, j-1);
            fval_2 = compute_special_feature(data, j-2);
            if(posLabel=='fwd')
                tempPred = quick_pred_fwd({fval_0, fval_1, fval_2});
            else
                tempPred = quick_pred_bkd({fval_0, fval_1, fval_2});
            end
            if(~isnan(tempPred))
                curPrediction(j) = tempPred;
                continue;
            end
            rawTesting = reshape(data(:, j-frameLen+1:j), [18*frameLen, 1]);
            curPrediction(j) = predict(svmModels, ((rawTesting'-trainingMean)./trainingStd)) + 0.5;
        end
        prediction = [prediction; curPrediction];
    end
end

function [labels, frameId, expId] = genLabels_test(alldata, testingIndxs, posLabel, startOffset, endOffset)
    assert(startOffset>=0 && endOffset>=0);
    labels = [];
    frameId = [];
    expId = [];
    switch(posLabel)
        case 'fwd'
            posLabel = 1;
        case 'bkd'
            posLabel = 3;
        otherwise
            posLabel = -1;
    end
    assert(posLabel~=-1);
    for i = 1:length(testingIndxs)
        curLabels = zeros(length(alldata(testingIndxs(i)).dFF_smooth), 1);
        tempFrameId = linspace(1, length(curLabels), length(curLabels));
        tempExpId = ones(size(curLabels))*testingIndxs(i);
        annos = alldata(testingIndxs(i)).anno;
        for j = 1:length(annos)
            startIndx = annos(j, 1)+1;
            endIndx = annos(j, 2)+1;
            curLabel = annos(j, 3);
            if(curLabel==posLabel)
                label = 1;
            elseif(curLabel==2 && posLabel==1)
                label = NaN;
            elseif(curLabel==4 && posLabel==3)
                label = NaN;
            elseif(curLabel==1 && posLabel==3)
                label = -1;
            elseif(curLabel==3 && posLabel==1)
                label = -1;
            else
                continue;
            end
            curLabels(startIndx:endIndx) = label;
            curLabels(startIndx:startIndx+startOffset) = NaN;
            curLabels(endIndx-endOffset:endIndx) = NaN;
        end
        labels = [labels; curLabels];
        frameId = [frameId; tempFrameId'];
        expId = [expId; tempExpId];
    end
end

%% functions copied

function fval = compute_special_feature(dff_smooth, nt)

    lt = 1:9;
    rt = 10:18;  % agreed upon convention (see runTCPserver.m)
    
    TRAIL = 20;
    ndx_stt = nt - TRAIL;
    ndx_end = nt;
    
    if ndx_stt<1; fval=0; return; end;
    
    [~,mxndx] = max(dff_smooth(:,ndx_stt:ndx_end),[],1);
    
    fval = sum(sign(diff(mxndx(lt)))) + sum(sign(diff(mxndx(rt))));
end

function result =  quick_pred_active(dFF_smooth, nt)
    result = 1;
    if max(abs(dFF_smooth(:,nt))) < 0.015
        result = 0.5;
    end
    if max(abs(dFF_smooth(:,nt))) < 0.01
        result = 0;
    end
end

function is_fwd = quick_pred_fwd(inprm)
% %     is_fwd = nan;
% %     if fval <= -5; is_fwd = 1; end
% %     if fval >= +5; is_fwd = 0; end
% %     
    fval = inprm{1};
    fval_o1 = inprm{2};
    fval_o2 = inprm{3};

    is_fwd = nan;
    if fval >= +5; is_fwd = 0; return; end
    if fval <= -7; 
        if fval_o1 < fval || fval_o2 < fval; return; end;
        if fval_o1 > fval && fval_o2 > fval_o1; is_fwd = 1; return ; end;    
    end
end

function is_bkw = quick_pred_bkd(inprm)
% %     is_bkw = nan;
% %     if fval >= +5; is_bkw = 1; end
% %     if fval <= -5; is_bkw = 0; end

    fval = inprm{1};
    fval_o1 = inprm{2};
    fval_o2 = inprm{3};

    is_bkw = nan;
    if fval <= -5; is_bkw = 0; return; end
    if fval >= +7; 
        if fval_o1 > fval || fval_o2 > fval; return; end;
        if fval_o1 < fval && fval_o2 < fval_o1; is_bkw = 1; return ; end;    
    end

end
