%% used to generate training data for each experiment

%% clean working space
% clear all;
% close all;
% clc;

%% load data
load alldata;

%% init
expNum = size(alldata, 2);
frameLen = 100;
posLabel = 'fwd';
startOffset = 0;
endOffset = 0;
channels = size(alldata(1).dFF_smooth, 1);

%% start loop
dataX = cell(expNum, 1);
dataY = cell(expNum, 1);
frameId = cell(expNum, 1);
fprintf(['processing starts! ' num2str(expNum) ' experiments to be processed.\n']);
for i = 1:expNum
    X = [];
    Y = [];
    fprintf(['now processing experiment ' num2str(i) ' ...\n']);
    [posIndx, negIndx] = genLabels(alldata, i, posLabel, startOffset, endOffset);
    posIndx = posIndx(posIndx > frameLen);
    negIndx = negIndx(negIndx > frameLen);
    [posIndx, negIndx] = downSample(1.0, 0.5, posIndx, negIndx);
    rawData = alldata(i).dFF_smooth;
    labels = [posIndx negIndx];
    dataY{i} = [ones(length(posIndx), 1); zeros(length(negIndx), 1)];
    frameId{i} = labels;
    for j = 1:length(labels)
        X = [X reshape(rawData(:, labels(j)-frameLen+1:labels(j)), [channels*frameLen, 1])];
    end
    dataX{i} = X';
end
fprintf('processing done!\n');

%% save results
save('dataX', 'dataX');
save('dataY', 'dataY');

%% for test


%% functions

% ------------------------------------------------------------------------------- %
% function 'genLabels' 
% ------------------------------------------------------------------------------- %
function [posIndx, negIndx] = genLabels(alldata, indx, posLabel, startOffset, endOffset)
    assert(startOffset>=0 && endOffset>=0);
    labels = ones(1, length(alldata(indx).dFF_smooth))*(-1);
    posIndx = [];
    negIndx = [];
    switch(posLabel)
        case 'fwd'
            posLabel = 1;
        case 'bkd'
            posLabel = 3;
        otherwise
            posLabel = -1;
    end
    assert(posLabel~=-1);
    for i = 1:length(alldata(indx).anno)
        startIndx = alldata(indx).anno(i, 1)+1;
        endIndx = alldata(indx).anno(i, 2)+1;
        curLabel = alldata(indx).anno(i, 3);
        if(curLabel==5)
            label = 0;
        else
            continue;
        end
        labels(startIndx:endIndx) = label;
    end
    for i = 1:length(alldata(indx).anno)
        startIndx = alldata(indx).anno(i, 1)+1;
        endIndx = alldata(indx).anno(i, 2)+1;
        curLabel = alldata(indx).anno(i, 3);
        if(curLabel==5)
            label = 0;
        else
            continue;
        end
        labels(startIndx:endIndx) = label;
    end
    for i = 1:length(alldata(indx).anno)
        startIndx = alldata(indx).anno(i, 1)+1;
        endIndx = alldata(indx).anno(i, 2)+1;
        curLabel = alldata(indx).anno(i, 3);
        if(curLabel==posLabel)
            label = 1;
        elseif(curLabel==2 && posLabel==1)
            label = -1;
        elseif(curLabel==4 && posLabel==3)
            label = -1;
        else
            continue;
        end
        labels(startIndx:endIndx) = label;
        labels(startIndx:startIndx+startOffset) = -1;
        labels(endIndx-endOffset:endIndx) = -1;
    end
    for i = 1:length(labels)
        if(labels(i)==1)
            posIndx = [posIndx i];
            continue;
        end
        if(labels(i)==0)
            negIndx = [negIndx i];
            continue;
        end
    end
end

% ------------------------------------------------------------------------------- %
% function 'downSample' 
% ------------------------------------------------------------------------------- %
function [posIndx, negIndx] = downSample(posKeepRate, negKeepRate, posIndx, negIndx)
    posIndx = posIndx(rand(size(posIndx)) <= posKeepRate);
    negIndx = negIndx(rand(size(negIndx)) <= negKeepRate);
end