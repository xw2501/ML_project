%% simple feature plot

%% load
load alldata;

%% 
orig_labels = [];
frame_id = [];
exp_id = [];
exp_num = length(alldata);
data_slc = cell(exp_num);
threshold = 0.015;

for i = 1:exp_num
    leng = size(alldata(i).dFF_smooth, 2);
    data = alldata(i).dFF_smooth;
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
    selected = find(max(abs(data), [], 1)>=threshold);
    labels = labels(selected);
    orig_labels = [orig_labels; labels];
    frame_id_tmp = (linspace(0, leng-1, leng))';
    frame_id_tmp = frame_id_tmp(selected);
    frame_id = [frame_id; frame_id_tmp];
    exp_id_tmp = ones(leng, 1) * i;
    exp_id_tmp = exp_id_tmp(selected);
    exp_id = [exp_id; exp_id_tmp];
    data_slc{i} = selected;
end

%% n_bound = 5, p_bound = 1:10
% len = 10;
% fwd_drop_fraction = zeros(len, 1);
% bkd_drop_fraction = zeros(len, 1);
% fwd_roc_pred = [];
% bkd_roc_pred = [];
% count = 1;
% for n_bound = 5:5
%     for p_bound = 1:10
%         fwd_selected = [];
%         bkd_selected = [];
%         for i = 1:exp_num
%             fwd_temp = nan(length(data_slc{i}), 1);
%             bkd_temp = nan(length(data_slc{i}), 1);
%             for j = 1:length(data_slc{i})
%                 data = alldata(i).dFF_smooth;
%                 fval_0 = compute_special_feature(data, data_slc{i}(j));
%                 fval_1 = compute_special_feature(data, data_slc{i}(j)-1);
%                 fval_2 = compute_special_feature(data, data_slc{i}(j)-2);
%                 if(quick_pred_fwd({fval_0, fval_1, fval_2}, n_bound, p_bound)==1)
%                     fwd_temp(j) = 1;
%                 elseif(quick_pred_fwd({fval_0, fval_1, fval_2}, n_bound, p_bound)==0)
%                     fwd_temp(j) = 0;
%                 else
%                     continue;
%                 end
%                 if(quick_pred_bkd({fval_0, fval_1, fval_2}, n_bound, p_bound)==1)
%                     bkd_temp(j) = 1;
%                 elseif(quick_pred_bkd({fval_0, fval_1, fval_2}, n_bound, p_bound)==0)
%                     bkd_temp(j) = 0;
%                 else
%                     continue;
%                 end
%             end
%             fwd_selected = [fwd_selected; fwd_temp];
%             bkd_selected = [bkd_selected; bkd_temp];
%         end
%         fwd_drop_fraction(count) = length([find(fwd_selected==1); find(fwd_selected==0)]) / length(fwd_selected);
%         bkd_drop_fraction(count) = length([find(bkd_selected==1); find(bkd_selected==0)]) / length(bkd_selected);
%         fwd_roc_pred = [fwd_roc_pred fwd_selected];
%         bkd_roc_pred = [bkd_roc_pred bkd_selected];
%         count = count + 1;
%     end
% end
% 
% figure, plot([1:10], fwd_drop_fraction, '-s');
% title('fwd nbound=5 pbound=1:10 drop fraction');
% savefig(gcf,  'nbound5_pbound1to10\fwd nbound5 pbound1to10 drop fraction.fig');
% figure, plot([1:10], bkd_drop_fraction, '-s');
% title('bkd nbound=5 pbound=1:10 drop fraction');
% savefig(gcf,  'nbound5_pbound1to10\bkd nbound5 pbound1to10 drop fraction.fig');
% pos_fwd = zeros(size(orig_labels));
% pos_fwd(find(orig_labels==5)) = 1;
% pos_bkd = zeros(size(orig_labels));
% pos_bkd(find(orig_labels==3)) = 1;
% 
% fwd_TP = zeros(10, 1);  bkd_TP = zeros(10, 1);
% fwd_FP = zeros(10, 1);  bkd_FP = zeros(10, 1);
% fwd_TN = zeros(10, 1);  bkd_TN = zeros(10, 1);
% fwd_FN = zeros(10, 1);  bkd_FN = zeros(10, 1);
% for i = 1:10
%     fwd_truth = pos_fwd(~isnan(fwd_roc_pred(:, i)));
%     fwd_pred = fwd_roc_pred(~isnan(fwd_roc_pred(:, i)), i);
%     bkd_truth = pos_bkd(~isnan(bkd_roc_pred(:, i)));
%     bkd_pred = bkd_roc_pred(~isnan(bkd_roc_pred(:, i)), i);
% 
%     fwd_TP(i) = length(find(fwd_truth(fwd_pred==1)==1));
%     fwd_FP(i) = length(find(fwd_truth(fwd_pred==1)==0));
%     fwd_TN(i) = length(find(fwd_truth(fwd_pred==0)==0));
%     fwd_FN(i) = length(find(fwd_truth(fwd_pred==0)==1));
%     
%     bkd_TP(i) = length(find(bkd_truth(bkd_pred==1)==1));
%     bkd_FP(i) = length(find(bkd_truth(bkd_pred==1)==0));
%     bkd_TN(i) = length(find(bkd_truth(bkd_pred==0)==0));
%     bkd_FN(i) = length(find(bkd_truth(bkd_pred==0)==1));
% end
% 
% fwd_TPR = fwd_TP./(fwd_TP+fwd_FN);
% fwd_FPR = fwd_FP./(fwd_TN+fwd_FP);
% bkd_TPR = bkd_TP./(bkd_TP+bkd_FN);
% bkd_FPR = bkd_FP./(bkd_TN+bkd_FP);
% 
% fwd_TPR(isnan(fwd_TPR)) = 0;
% fwd_FPR(isnan(fwd_FPR)) = 0;
% bkd_TPR(isnan(bkd_TPR)) = 0;
% bkd_FPR(isnan(bkd_FPR)) = 0;
% 
% figure, plot(fwd_FPR, fwd_TPR, '-s');
% xlabel('FPR');
% ylabel('TPR');
% title('ROC of forward classifier');
% axis([0 1 0 1]);
% text(fwd_FPR(1:2:7)+0.03, fwd_TPR(1:2:7)-0.03, {'1', '3', '5', '7'});
% text(fwd_FPR(9)+0.03, fwd_TPR(9)+0.03, '9');
% savefig(gcf,  'nbound5_pbound1to10\forward ROC nbound5 pbound1to10.fig');
% 
% figure, plot(bkd_FPR, bkd_TPR, '-s');
% xlabel('FPR');
% ylabel('TPR');
% title('ROC of backward classifier');
% axis([0 1 0 1]);
% text(bkd_FPR(1:2:7), bkd_TPR(1:2:7)-0.03, {'1', '3', '5', '7'});
% text(bkd_FPR(9), bkd_TPR(9)-0.03, '9');
% savefig(gcf,  'nbound5_pbound1to10\backward ROC nbound5 pbound1to10.fig');
% 
% fwd_fwd_fraction = zeros(len, 1);
% fwd_bkd_fraction = zeros(len, 1);
% fwd_act_fraction = zeros(len, 1);
% fwd_sln_fraction = zeros(len, 1);
% bkd_fwd_fraction = zeros(len, 1);
% bkd_bkd_fraction = zeros(len, 1);
% bkd_act_fraction = zeros(len, 1);
% bkd_sln_fraction = zeros(len, 1);
% for i = 1:len
%     fwd_fwd_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==5)) / length(find(isnan(fwd_roc_pred(:, i))==0));
%     fwd_bkd_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==3)) / length(find(isnan(fwd_roc_pred(:, i))==0));
%     fwd_act_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==1)) / length(find(isnan(fwd_roc_pred(:, i))==0));
%     fwd_sln_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==0)) / length(find(isnan(fwd_roc_pred(:, i))==0));
%     bkd_fwd_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==5)) / length(find(isnan(bkd_roc_pred(:, i))==0));
%     bkd_bkd_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==3)) / length(find(isnan(bkd_roc_pred(:, i))==0));
%     bkd_act_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==1)) / length(find(isnan(bkd_roc_pred(:, i))==0));
%     bkd_sln_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==0)) / length(find(isnan(bkd_roc_pred(:, i))==0));
% end
% figure,
% hold on;
% plot([1:len], fwd_fwd_fraction, '-s');
% plot([1:len], fwd_fwd_fraction+fwd_bkd_fraction, '-o');
% plot([1:len], fwd_fwd_fraction+fwd_bkd_fraction+fwd_act_fraction, '-*');
% plot([1:len], fwd_fwd_fraction+fwd_bkd_fraction+fwd_act_fraction+fwd_sln_fraction, '-d');
% title('fwd classifier cumulate data fraction');
% legend('fwd', 'fwd+bkd', 'fwd+bkd+act', 'one');
% savefig(gcf, 'nbound5_pbound1to10\fwd classifier cumulate data fraction');
% figure,
% hold on;
% plot([1:len], bkd_fwd_fraction, '-s');
% plot([1:len], bkd_fwd_fraction+bkd_bkd_fraction, '-o');
% plot([1:len], bkd_fwd_fraction+bkd_bkd_fraction+bkd_act_fraction, '-*');
% plot([1:len], bkd_fwd_fraction+bkd_bkd_fraction+bkd_act_fraction+bkd_sln_fraction, '-d');
% title('bkd classifier cumulate data fraction');
% legend('fwd', 'fwd+bkd', 'fwd+bkd+act', 'one');
% savefig(gcf, 'nbound5_pbound1to10\bkd classifier cumulate data fraction');

%% n_bound = 1:10, p_bound = 6
len = 10;
fwd_drop_fraction = zeros(len, 1);
bkd_drop_fraction = zeros(len, 1);
fwd_roc_pred = [];
bkd_roc_pred = [];
count = 1;
for p_bound = 6:6
    for n_bound = 1:10
        fwd_selected = [];
        bkd_selected = [];
        for i = 1:exp_num
            fwd_temp = nan(length(data_slc{i}), 1);
            bkd_temp = nan(length(data_slc{i}), 1);
            for j = 1:length(data_slc{i})
                data = alldata(i).dFF_smooth;
                fval_0 = compute_special_feature(data, data_slc{i}(j));
                fval_1 = compute_special_feature(data, data_slc{i}(j)-1);
                fval_2 = compute_special_feature(data, data_slc{i}(j)-2);
                if(quick_pred_fwd({fval_0, fval_1, fval_2}, n_bound, p_bound)==1)
                    fwd_temp(j) = 1;
                elseif(quick_pred_fwd({fval_0, fval_1, fval_2}, n_bound, p_bound)==0)
                    fwd_temp(j) = 0;
                else
                    continue;
                end
                if(quick_pred_bkd({fval_0, fval_1, fval_2}, n_bound, p_bound)==1)
                    bkd_temp(j) = 1;
                elseif(quick_pred_bkd({fval_0, fval_1, fval_2}, n_bound, p_bound)==0)
                    bkd_temp(j) = 0;
                else
                    continue;
                end
            end
            fwd_selected = [fwd_selected; fwd_temp];
            bkd_selected = [bkd_selected; bkd_temp];
        end
        fwd_drop_fraction(count) = length([find(fwd_selected==1); find(fwd_selected==0)]) / length(fwd_selected);
        bkd_drop_fraction(count) = length([find(bkd_selected==1); find(bkd_selected==0)]) / length(bkd_selected);
        fwd_roc_pred = [fwd_roc_pred fwd_selected];
        bkd_roc_pred = [bkd_roc_pred bkd_selected];
        count = count + 1;
    end
end

figure, plot([1:10], fwd_drop_fraction, '-s');
title('fwd pbound=6 nbound=1:10 drop fraction');
savefig(gcf,  'nbound1to10_pbound6\fwd pbound6 nbound1to10 drop fraction.fig');
figure, plot([1:10], bkd_drop_fraction, '-s');
title('bkd pbound=6 nbound=1:10 drop fraction');
savefig(gcf,  'nbound1to10_pbound6\bkd pbound6 nbound1to10 drop fraction.fig');
pos_fwd = zeros(size(orig_labels));
pos_fwd(find(orig_labels==5)) = 1;
pos_bkd = zeros(size(orig_labels));
pos_bkd(find(orig_labels==3)) = 1;

fwd_TP = zeros(10, 1);  bkd_TP = zeros(10, 1);
fwd_FP = zeros(10, 1);  bkd_FP = zeros(10, 1);
fwd_TN = zeros(10, 1);  bkd_TN = zeros(10, 1);
fwd_FN = zeros(10, 1);  bkd_FN = zeros(10, 1);
for i = 1:10
    fwd_truth = pos_fwd(~isnan(fwd_roc_pred(:, i)));
    fwd_pred = fwd_roc_pred(~isnan(fwd_roc_pred(:, i)), i);
    bkd_truth = pos_bkd(~isnan(bkd_roc_pred(:, i)));
    bkd_pred = bkd_roc_pred(~isnan(bkd_roc_pred(:, i)), i);

    fwd_TP(i) = length(find(fwd_truth(fwd_pred==1)==1));
    fwd_FP(i) = length(find(fwd_truth(fwd_pred==1)==0));
    fwd_TN(i) = length(find(fwd_truth(fwd_pred==0)==0));
    fwd_FN(i) = length(find(fwd_truth(fwd_pred==0)==1));
    
    bkd_TP(i) = length(find(bkd_truth(bkd_pred==1)==1));
    bkd_FP(i) = length(find(bkd_truth(bkd_pred==1)==0));
    bkd_TN(i) = length(find(bkd_truth(bkd_pred==0)==0));
    bkd_FN(i) = length(find(bkd_truth(bkd_pred==0)==1));
end

fwd_TPR = fwd_TP./(fwd_TP+fwd_FN);
fwd_FPR = fwd_FP./(fwd_TN+fwd_FP);
bkd_TPR = bkd_TP./(bkd_TP+bkd_FN);
bkd_FPR = bkd_FP./(bkd_TN+bkd_FP);

fwd_TPR(isnan(fwd_TPR)) = 0;
fwd_FPR(isnan(fwd_FPR)) = 0;
bkd_TPR(isnan(bkd_TPR)) = 0;
bkd_FPR(isnan(bkd_FPR)) = 0;

figure, plot(fwd_FPR, fwd_TPR, '-s');
xlabel('FPR');
ylabel('TPR');
title('ROC of forward classifier');
axis([0 1 0 1]);
text(fwd_FPR(1:2:7)+0.03, fwd_TPR(1:2:7)-0.03, {'1', '3', '5', '7'});
text(fwd_FPR(9)+0.03, fwd_TPR(9)+0.03, '9');
savefig(gcf,  'nbound1to10_pbound6\forward ROC pbound6 nbound1to10.fig');

figure, plot(bkd_FPR, bkd_TPR, '-s');
xlabel('FPR');
ylabel('TPR');
title('ROC of backward classifier');
axis([0 1 0 1]);
text(bkd_FPR(1:2:7), bkd_TPR(1:2:7)-0.03, {'1', '3', '5', '7'});
text(bkd_FPR(9), bkd_TPR(9)-0.03, '9');
savefig(gcf,  'nbound1to10_pbound6\backward ROC pbound6 nbound1to10.fig');

fwd_fwd_fraction = zeros(len, 1);
fwd_bkd_fraction = zeros(len, 1);
fwd_act_fraction = zeros(len, 1);
fwd_sln_fraction = zeros(len, 1);
bkd_fwd_fraction = zeros(len, 1);
bkd_bkd_fraction = zeros(len, 1);
bkd_act_fraction = zeros(len, 1);
bkd_sln_fraction = zeros(len, 1);
for i = 1:len
    fwd_fwd_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==5)) / length(find(isnan(fwd_roc_pred(:, i))==0));
    fwd_bkd_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==3)) / length(find(isnan(fwd_roc_pred(:, i))==0));
    fwd_act_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==1)) / length(find(isnan(fwd_roc_pred(:, i))==0));
    fwd_sln_fraction(i) = length(find(orig_labels(~isnan(fwd_roc_pred(:, i)))==0)) / length(find(isnan(fwd_roc_pred(:, i))==0));
    bkd_fwd_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==5)) / length(find(isnan(bkd_roc_pred(:, i))==0));
    bkd_bkd_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==3)) / length(find(isnan(bkd_roc_pred(:, i))==0));
    bkd_act_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==1)) / length(find(isnan(bkd_roc_pred(:, i))==0));
    bkd_sln_fraction(i) = length(find(orig_labels(~isnan(bkd_roc_pred(:, i)))==0)) / length(find(isnan(bkd_roc_pred(:, i))==0));
end
figure,
hold on;
plot([1:len], fwd_fwd_fraction, '-s');
plot([1:len], fwd_fwd_fraction+fwd_bkd_fraction, '-o');
plot([1:len], fwd_fwd_fraction+fwd_bkd_fraction+fwd_act_fraction, '-*');
plot([1:len], fwd_fwd_fraction+fwd_bkd_fraction+fwd_act_fraction+fwd_sln_fraction, '-d');
title('fwd classifier cumulate data fraction');
legend('fwd', 'fwd+bkd', 'fwd+bkd+act', 'one');
savefig(gcf, 'nbound1to10_pbound6\fwd classifier cumulate data fraction');
figure,
hold on;
plot([1:len], bkd_fwd_fraction, '-s');
plot([1:len], bkd_fwd_fraction+bkd_bkd_fraction, '-o');
plot([1:len], bkd_fwd_fraction+bkd_bkd_fraction+bkd_act_fraction, '-*');
plot([1:len], bkd_fwd_fraction+bkd_bkd_fraction+bkd_act_fraction+bkd_sln_fraction, '-d');
title('bkd classifier cumulate data fraction');
legend('fwd', 'fwd+bkd', 'fwd+bkd+act', 'one');
savefig(gcf, 'nbound1to10_pbound6\bkd classifier cumulate data fraction');

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