

%% load
%numrois = 6;
numrois = 18;

base_dir = '/groups/branson/home/verman/data/work/projs/chen';
data_dir = sprintf('%s/data/%0dROIs',base_dir,numrois);

alldata = class_process_data.load_and_prep_data(base_dir, data_dir, numrois);

%% trace normalization

ndatasets = length(alldata);
for dsi = 1:ndatasets
fprintf('normalizing dataset %d/%d\n', dsi,ndatasets);
  tm = alldata(dsi).trace_mean';
  [T,nroi] = size(tm);
  nm = nan(size(tm));
  nms = nan(size(tm));


  x = 700;
  p = 5;
  e = 30;
  b = 100;

  filtersd = 4;
  for t=1:T
      [nm(t,:), nms(t,:)] = preprocess_trace('F',tm, 'dFF', nm, 'timeIndex', t, ...
                                             'windowsize', x, 'percentile', p, 'epsilon', e, 'background', b, ...
                                             'filtersd',filtersd);
  end
 
  alldata(dsi).dFF = nm';
  alldata(dsi).dFF_smooth = nms';
 
end



return;


%% visualize -- each dataset max roi

alldata;
datasource = 'dFF_smooth';

lowthresh = -1; %0.015;


figure; 
hold on;
for dsi = 1:length(alldata)
    ds = alldata(dsi);
    t = max(abs(ds.(datasource)),[],1);
    t(t<lowthresh) = nan;
    
    fprintf('%d: %1.2f\n',dsi,sum(isnan(t))/length(t));
    
    n = ds.name;
    plot(t,'-');
    pause;
end
title('max roi','interpreter','none');


%% visualize -- each dataset across rois

alldata;
datasource = 'dFF_smooth';


for dsi = 1:length(alldata)
    ds = alldata(dsi);
    t = ds.(datasource);
    n = ds.name;
    
    figure; 
    hold on;
    
    for roii=1:numrois
        plot(t(roii,:),'-');
    end  
    title(sprintf('dataset: %s',n),'interpreter','none');
    pause;
end

%% visualize -- each roi across datasets
alldata;
datasource = 'dFF_smooth';

for roii=1:numrois
    figure; 
    hold on;
    
    for dsi = 1:length(alldata)
        ds = alldata(dsi);
        
        t = ds.(datasource);
        plot(t(roii,:),'-');
    end  
    title(sprintf('roi: %d',roii));
    pause;
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%   OLD STUFF  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%i1 = 1;
%i2 = 1;

for j=1%:9*2
    i1 = j; %[i1,i2] = ind2sub([9 2],j);
figure;
hold on;
lstr = '';
%tstr = sprintf('A%d, D%d',i1,i2);
tstr = sprintf('J = %d',j);
for i=1:length(alldata)
    ad = alldata(i);
    %t = squeeze(ad.trace_mean(i1,i2,:));
    %n = squeeze(ad.trace_normalized(i1,i2,:));
    t = ad.trace_mean(j,:);
    n = ad.trace_normalized(j,:);
    plot(n);
    lstr = sprintf('%s, ''%s''',lstr,ad.name);
end


eval(sprintf('legend(%s);',strrep(lstr(2:end),'_',' ')));
xlabel('time (in frames)');
ylabel('calcium activity level');
title(tstr);
xl = xlim();
plot(xl,[0 0],'k-','linewidth',2);
end

%%
% % 
% % %% create normalized data
% % for triali=1:length(alldata)
% %     ad = alldata(triali);
% %     dsz=size(ad.trace_mean);
% %     trace_normalized = nan(dsz);
% %     trace_mean = ad.trace_mean;
% %     for i=1:dsz(1)*dsz(2)
% %         [a1,a2]=ind2sub(dsz(1:2),i);
% %         t = squeeze(trace_mean(a1,a2,:));
% %         b=nan(size(t));
% %         for j=1:length(t)
% %             b(j) = min(t(1:j));
% %         end
% %         n = t-b;
% %         trace_normalized(a1,a2,:) = n;
% %     end
% %     alldata(triali).trace_normalized = trace_normalized;
% % end

%%  extract training data (positive label ==  forward)
poslabel_name = 'fwd';
tslen = 100;
%datasource = 'trame_mean';
datasource = 'trace_normalized';

posstep=1;
negstep=1;

[all_X,all_Y,trial_id,frame_id, traialnames] = ...
               class_process_data.extract_training_data(alldata,poslabel_name,'tslen',tslen, 'datasource',datasource, ...
                                                        'posstep',posstep,'negstep',negstep ...
               );

%% train/test split
train_frac = 0.5;
[train_X, train_Y, test_X,test_Y, trainidx] = class_process_data.train_test_split(all_X,all_Y, trial_id, frame_id,train_frac);

%%  SVM train data

[normdat,mn,st] = standardizeCols(train_X);
odat = normdat;
labels = train_Y;
    
svm_cv_prm = {'train_test_balance',[true true], 'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto'};
fwd_svm_results = [];
[fwd_svm_results.wt_vec,fwd_svm_results.svm_struct,fwd_svm_results.TRloss,fwd_svm_results.predy] = class_svm.svm(odat,labels, svm_cv_prm);
fwd_svm_results.train_frame_id = frame_id(trainidx);
fwd_svm_results.train_trial_id = trial_id(trainidx);
fwd_svm_results.mn = mn;
fwd_svm_results.st = st;

%%  SVM cross validation
trial_id;

[normdat,mn,st] = standardizeCols(x);
odat = normdat;
labels = y;
svm_cv_prm = {'train_test_balance',[true true], 'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto'};


nfolds = 3;
npts = length(trial_id);

cvndx = cross_validate_split_bygroup(npts,nfolds,trial_id);
[min_avgTSloss, all_TSloss, min_avgTRloss, all_TRloss] = class_svm.svm_cv(odat, labels, svm_cv_prm, cvndx)   


%% full prediction
datasource;

for i=1:length(alldata)
    ad = alldata(i);
    pred = [];
    [pred.x,pred.frame_id,pred.f_desc] = class_process_data.all_instance_data(ad.(datasource),tslen);
    nx = size(pred.x,1);
    pred.norm_x = (pred.x - repmat(fwd_svm_results.mn,nx,1)) ./ repmat(fwd_svm_results.st,nx,1);
    pred.plabel = class_svm.pred_fxn(pred.norm_x, struct('svmstruct',fwd_svm_results.svm_struct));
    alldata(i).pred = pred;
end



%% visualize prediction



alldata;
svm_results;
labels;

posid = class_process_data.resolvelbl(poslabel_name);

i1 = 1;
i2 = 1;


for i=1:length(alldata)
    figure; 
    hold on;
    
    ad = alldata(i);
    t = squeeze(ad.trace_mean(i1,i2,:));
    
    tidx = svm_results.train_trial_id==i;
    %pi = svm_results.predy(tidx);
    trainfi = svm_results.train_frame_id(tidx);
    allfi = ad.pred.frame_id;
    allpi = ad.pred.plabel;
    
    pani = ad.anno(ad.anno(:,end)==posid,:);
    
    plot(t);
    if ~isempty(pani)
        for j=1:size(pani,1)
            plot([pani(j,1) pani(j,2)], [0 0], '-','linewidth', 2,'color',[0 0 0]);
        end
    end
    plot(trainfi,zeros(size(trainfi)),'ok');
    plot(allfi(allpi==0),zeros(size(allfi(allpi==0))),'xb');
    plot(allfi(allpi==1),zeros(size(allfi(allpi==1))),'xr');
    
    title(ad.name,'interpreter','none');
end


%%  SVM all data

[normdat,mn,st] = standardizeCols(all_X);
odat = normdat;
labels = all_Y;
    
svm_cv_prm = {'train_test_balance',[true true], 'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto'};
fwd_svm_results = [];
[fwd_svm_results.wt_vec,fwd_svm_results.svm_struct,fwd_svm_results.TRloss,fwd_svm_results.predy] = class_svm.svm(odat,labels, svm_cv_prm);
%fwd_svm_results.train_frame_id = frame_id(trainidx);
%fwd_svm_results.train_trial_id = trial_id(trainidx);
fwd_svm_results.mn = mn;
fwd_svm_results.st = st;

%% save fwd weight params
savefilename = 'wtparams_fwd.mat';
wtprms_fwdwave = [];
wtprms_fwdwave_start = [];

wtprms_fwdwave.tslen = tslen;
wtprms_fwdwave.svm_prms = fwd_svm_results;



save(savefilename, 'wtprms_fwdwave', 'wtprms_fwdwave_start');





%%




