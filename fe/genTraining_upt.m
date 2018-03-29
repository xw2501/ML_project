%% load data
load alldata;

%% init
expNum = size(alldata, 2);
frameLen = 100;
posLabel = 1;   % 1-forward, 3-backward
startOffset = 0;
endOffset = 0;
threshold = 0.015;
n_bound = 5;
p_bound = 6;
channels = size(alldata(1).dFF_smooth, 1);

%% start loop

dataX = cell(expNum, 1);
dataY = cell(expNum, 1);
frameId = cell(expNum, 1);

for i = 1:expNum
    data = alldata(i).dFF_smooth;
    labels = zeros(length(data), 1);
    
    % remove silence state
    labels(max(abs(data))>=threshold) = 1;
    selected = find(labels==1);
    
    % remove data through feature
    for j = 1:length(selected)
        fval_0 = compute_special_feature(data, selected(i));
        fval_1 = compute_special_feature(data, selected(i)-1);
        fval_2 = compute_special_feature(data, selected(i)-2);
        
        if(posLabel==1)
            res = quick_pred_fwd({fval_0,fval_1,fval_2}, n_bound, p_bound);
            if(~isnan(res))
                labels(selected(i)) = 0;
            end
        elseif(posLabel==3)
            res = quick_pred_bkd({fval_0,fval_1,fval_2}, n_bound, p_bound);
            if(~isnan(res))
                labels(selected(i)) = 0;
            end
        else
            error('wrong label!');
        end
    end
    selected = find(labels==1);
    selected = selected(find(selected>frameLen));
    
    % make labels
    labels = ones(size(labels)) * (-1);
    annos = alldata(i).anno;
    for j = length(annos):-1:1
        label = annos(j, 3);
        startIndx = annos(j, 1)+1;
        endIndx = annos(j, 2)+1;
        if(label==6)                % skip silence state
            continue;
        elseif(label==5)            % label 0 for active
            label = 0;
            labels(startIndx:endIndx) = label;
        elseif(label==posLabel)     % label 1 for positive
            label = 1;
            labels(startIndx:endIndx) = label;
            labels(startIndx:startIndx+startOffset) = -1;
            labels(endIndx-endOffset:endIndx) = -1;
        elseif(label==(posLabel+1))   % unlabel weak states
            label = -1;
            labels(startIndx:endIndx) = label;
        else                        % label 0 for the rest states
            label = 0;
            labels(startIndx:endIndx) = label;
            labels(startIndx:startIndx+startOffset) = -1;
            labels(endIndx-endOffset:endIndx) = -1;
        end
    end
    
    % downsample
    posIndx = selected(find(labels(selected)==1));
    negIndx = selected(find(labels(selected)==0));
    [posIndx, negIndx] = downSample(1.0, 0.5, posIndx, negIndx);
    selected = [posIndx; negIndx];
    dataY{i} = [ones(length(posIndx), 1); zeros(length(negIndx), 1)];
    X = [];
    for j = 1:length(selected)
        X = [X reshape(data(:, selected(j)-frameLen+1:selected(j)), [channels*frameLen, 1])];
    end
    dataX{i} = X';
    frameId{i} = selected;
end

save('dataX', 'dataX');
save('dataY', 'dataY');
save('frameId', 'frameId');

%% functions

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

function is_fwd = quick_pred_fwd(inprm, n_bound, p_bound)
% %     is_fwd = nan;
% %     if fval <= -5; is_fwd = 1; end
% %     if fval >= +5; is_fwd = 0; end
% %     
    fval = inprm{1};
    fval_o1 = inprm{2};
    fval_o2 = inprm{3};

    is_fwd = nan;
    if fval >= n_bound; is_fwd = 0; return; end
    if fval <= (-1*p_bound); 
        if fval_o1 < fval || fval_o2 < fval; return; end;
        if fval_o1 > fval && fval_o2 > fval_o1; is_fwd = 1; return ; end;    
    end
end

function is_bkw = quick_pred_bkd(inprm, n_bound, p_bound)
% %     is_bkw = nan;
% %     if fval >= +5; is_bkw = 1; end
% %     if fval <= -5; is_bkw = 0; end

    fval = inprm{1};
    fval_o1 = inprm{2};
    fval_o2 = inprm{3};

    is_bkw = nan;
    if fval <= (-1*n_bound); is_bkw = 0; return; end
    if fval >= p_bound; 
        if fval_o1 > fval || fval_o2 > fval; return; end;
        if fval_o1 < fval && fval_o2 < fval_o1; is_bkw = 1; return ; end;    
    end

end

function [posIndx, negIndx] = downSample(posKeepRate, negKeepRate, posIndx, negIndx)
    posIndx = posIndx(rand(size(posIndx)) <= posKeepRate);
    negIndx = negIndx(rand(size(negIndx)) <= negKeepRate);
end