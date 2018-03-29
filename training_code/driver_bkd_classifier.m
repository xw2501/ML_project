alldata;



%% extract bkw wave ndxs;

% display_result = true;
display_result = false;

plbls = {'bkd'};
pin = {'st0' 'en0'};
pig = {'st-5' 'en+5'};
nlbls = {}; 
nin = {}; 
nig = {};
ilbls = {'weak_bkd'};
iin = {'st0' 'en0'};

[pos_ndxs,neg_ndxs,ign_ndxs] = class_process_data.get_ndxs(alldata,plbls,pin,pig, nlbls,nin,nig, ilbls,iin, display_result);



%% create training data
tslen = 100;
datasource = 'dFF_smooth';
point_density_pos = 1;
point_density_neg = .5;

min_intensity_thresh = 0.015;

[x,y,trial_id,frame_id, traialnames, f_desc] = class_process_data.extract_training_generic(alldata,pos_ndxs,neg_ndxs, ...
                                                     'tslen',tslen, 'datasource', datasource, 'point_density_pos', point_density_pos, 'point_density_neg', point_density_neg, ...
                                                      'min_intensity_thresh', min_intensity_thresh);


%%

%%  SVM all data

[normdat,mn,st] = standardizeCols(x);
%[pca_proj_matrix,odat, ~,~,pca_var_explained] = pca(normdat, 'numcomponents',100);
odat = normdat; pca_proj_matrix = [];
labels = y;
%%    
% % svm_cv_prm = {'train_test_balance',[true true], 'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto'};
% % bkd_svm_results = [];
% % [bkd_svm_results.wt_vec,bkd_svm_results.svm_struct,bkd_svm_results.TRloss,bkd_svm_results.predy] = class_svm.svm(odat,labels, svm_cv_prm);
% % bkd_svm_results.train_frame_id = frame_id; %(trainidx);
% % bkd_svm_results.train_trial_id = trial_id; %(trainidx);
% % bkd_svm_results.min_intensity_thresh = min_intensity_thresh;
% % bkd_svm_results.mn = mn;
% % bkd_svm_results.st = st;
% % bkd_svm_results.pca_proj_matrix = pca_proj_matrix;
% % 
% % bkd_svm_results.TRloss

%%    
svm_cv_prm = {'train_test_balance',[true true], 'Learner','svm', 'lambda', 0.5, 'cost', [0 10; 1 0]};
bkd_svm_results = [];
[bkd_svm_results.wt_vec,bkd_svm_results.svm_struct,bkd_svm_results.TRloss,bkd_svm_results.predy] = class_linear_predictor.linpredict(odat,labels, svm_cv_prm);
bkd_svm_results.train_frame_id = frame_id; %(trainidx);
bkd_svm_results.train_trial_id = trial_id; %(trainidx);
bkd_svm_results.min_intensity_thresh = min_intensity_thresh;
bkd_svm_results.mn = mn;
bkd_svm_results.st = st;
bkd_svm_results.pca_proj_matrix = pca_proj_matrix;

bkd_svm_results.TRloss
%% save active weight params
savefilename = 'wtparams_bkw.mat';
wtprms_bkdwave = [];

wtprms_bkdwave.tslen = tslen;
wtprms_bkdwave.svm_prms = bkd_svm_results;



save(savefilename, 'wtprms_bkdwave');


