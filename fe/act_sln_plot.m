%% 'act'-'sln' threshold
% note:
%   in experiment 24, state labels from index 1479 to 1578 missing

%% load
% load alldata;

%% 
orig_labels = [];
frame_id = [];
exp_id = [];
exp_num = length(alldata);
posLabel = 0;       % 0 for silence, 1 for activity

%% origin label
for i = 1:exp_num
    leng = size(alldata(i).dFF_smooth, 2);
    labels = nan(leng, 1);
    for j = 1:length(alldata(i).anno)
        label = alldata(i).anno(j, 3);
        if(label~=5 && label~=6)
            continue;
        end
        if(label==posLabel)
            label = 6 - label;
        else
            label = label - 5;
        end
        startIndx = alldata(i).anno(j, 1)+1;
        endIndx = alldata(i).anno(j, 2)+1;
        labels(startIndx:endIndx) = label;
    end
    orig_labels = [orig_labels; labels];
    frame_id = [frame_id; (linspace(0, leng-1, leng))'];
    exp_id = [exp_id; ones(leng, 1) * i];
end

%% threshold
threshold_start = 0.0005;
threshold_end = 0.05;
step = 0.0005;
threshold = threshold_start:step:threshold_end;
TP = zeros(length(threshold), 1);
FP = zeros(length(threshold), 1);
TN = zeros(length(threshold), 1);
FN = zeros(length(threshold), 1);
pred_labels = [];

for j = 1:length(threshold)
    temp_labels = [];
    for i = 1:exp_num
        data = alldata(i).dFF_smooth;
        labels = zeros(size(data, 2), 1);
        if(posLabel==1)
            labels(max(abs(data))>=threshold(j)) = 1;
        else
            labels(max(abs(data))<threshold(j)) = 1;
        end
        temp_labels = [temp_labels; labels];
    end
    pred_labels = [pred_labels temp_labels];
    TP(j) = length(find(orig_labels(temp_labels==1)==1));    % turth 1, pred 1
    FP(j) = length(find(orig_labels(temp_labels==1)==0));    % turth 0, pred 1
    TN(j) = length(find(orig_labels(temp_labels==0)==0));    % truth 0, pred 0
    FN(j) = length(find(orig_labels(temp_labels==0)==1));    % truth 1, pred 0
end

TPR = TP./(TP+FN);
FPR = FP./(TN+FP);
figure, plot(FPR, TPR, '-bs');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC of threshold');
axis([0 1 0 1]);
hold on;
text(0.5, 0.7, 'threshold starting from 0.0005');
text(0.5, 0.65, 'at left bottom corner');
text(0.5, 0.6, 'threshold ending at 0.05');
text(0.5, 0.55, 'threshold step size is 0.0005');

%%
drop_fraction = zeros(length(threshold), 1);
fwd_drop = zeros(size(drop_fraction));
bkd_drop = zeros(size(drop_fraction));
act_drop = zeros(size(drop_fraction));
sln_drop = zeros(size(drop_fraction));
fwd_fraction = zeros(size(drop_fraction));
bkd_fraction = zeros(size(drop_fraction));
act_fraction = zeros(size(drop_fraction));
sln_fraction = zeros(size(drop_fraction));

orig_labels = [];
for i = 1:exp_num
    leng = size(alldata(i).dFF_smooth, 2);
    labels = nan(leng, 1);
    for j = length(alldata(i).anno):-1:1
        label = alldata(i).anno(j, 3);
        if(label==2 || label==4)
            label = 5;
        end
        % 'fwd' label = 5
        % 'bkd' label = 3
        % 'act' label = 1
        % 'sln' label = 0
        label = 6 - label;
        startIndx = alldata(i).anno(j, 1)+1;
        endIndx = alldata(i).anno(j, 2)+1;
        labels(startIndx:endIndx) = label;
    end
    orig_labels = [orig_labels; labels];
end

for i = 1:length(threshold)
    selected = find(pred_labels(:, i)==0);
    drop_fraction(i) = (length(orig_labels)-length(selected)) / length(orig_labels);
    fwd_drop(i) = length(find(orig_labels(setdiff(linspace(1, length(orig_labels), length(orig_labels)), selected))==5)) / length(find(orig_labels==5));
    bkd_drop(i) = length(find(orig_labels(setdiff(linspace(1, length(orig_labels), length(orig_labels)), selected))==3)) / length(find(orig_labels==3));
    act_drop(i) = length(find(orig_labels(setdiff(linspace(1, length(orig_labels), length(orig_labels)), selected))==1)) / length(find(orig_labels==1));
    sln_drop(i) = length(find(orig_labels(setdiff(linspace(1, length(orig_labels), length(orig_labels)), selected))==0)) / length(find(orig_labels==0));
    fwd_fraction(i) = length(find(orig_labels(selected)==5)) / length(selected);
    bkd_fraction(i) = length(find(orig_labels(selected)==3)) / length(selected);
    act_fraction(i) = length(find(orig_labels(selected)==1)) / length(selected);
    sln_fraction(i) = length(find(orig_labels(selected)==0)) / length(selected);
end

figure, plot(threshold, drop_fraction, '-s');
hold on;
plot(threshold, fwd_drop, '-o');
plot(threshold, bkd_drop, '-*');
plot(threshold, act_drop, '-d');
plot(threshold, sln_drop, '-^');
title('data drop fraction over threshold');
legend('alldata', 'forward', 'backward', 'active', 'silence');
xlabel('threshold');
ylabel('fraction');

figure,
hold on;
plot(threshold, fwd_fraction, '-s');
plot(threshold, fwd_fraction+bkd_fraction, '-o');
plot(threshold, fwd_fraction+bkd_fraction+act_fraction, '-*');
plot(threshold, ones(size(threshold)), '-d');
title('cumulate data fraction');
legend('fwd', 'fwd+bkd', 'fwd+bkd+act', 'one');
xlabel('threshold');
ylabel('fraction');

figure,
hold on;
plot(threshold, fwd_fraction, '-s');
plot(threshold, bkd_fraction, '-o');
plot(threshold, act_fraction, '-*');
plot(threshold, sln_fraction, '-d');
title('individual data fraction');
legend('fwd', 'bkd', 'act', 'sln');
xlabel('threshold');
ylabel('fraction');
