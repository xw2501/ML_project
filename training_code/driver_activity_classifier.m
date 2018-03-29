%%
alldata;



%% extract activity ndxs;

display_result = true;
%display_result = false;

plbls = {'activity'};
pin = {'st0' 'en0'};
pig = {'st-5' 'en+5'};
nlbls = {}; 
nin = {}; 
nig = {};
ilbls = {}; %'weak_fwd'};
iin = {}; %'st0' 'en0'};

[pos_ndxs,neg_ndxs,ign_ndxs] = class_process_data.get_ndxs(alldata,plbls,pin,pig, nlbls,nin,nig, ilbls,iin, display_result);



%% create training data

tslen = 100;
datasource = 'dFF_smooth';
point_density = 1;
min_intensity_thresh = 0.015;

[x,y,trial_id,frame_id, traialnames, f_desc] = class_process_data.extract_training_generic(alldata,pos_ndxs,neg_ndxs, ...
                                                     'tslen',tslen, 'datasource', datasource, 'point_density', point_density, ...
                                                      'min_intensity_thresh', min_intensity_thresh);


%%

%%  SVM all data

[normdat,mn,st] = standardizeCols(x);
%[pca_proj_matrix,odat, ~,~,pca_var_explained] = pca(normdat, 'numcomponents',100);
odat = normdat; pca_proj_matrix = [];
labels = y;
%%    
svm_cv_prm = {'train_test_balance',[true true], 'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto'};
act_svm_results = [];
[act_svm_results.wt_vec,act_svm_results.svm_struct,act_svm_results.TRloss,act_svm_results.predy] = class_svm.svm(odat,labels, svm_cv_prm);
act_svm_results.train_frame_id = frame_id; %(trainidx);
act_svm_results.train_trial_id = trial_id; %(trainidx);
act_svm_results.min_intensity_thresh = min_intensity_thresh;
act_svm_results.mn = mn;
act_svm_results.st = st;
act_svm_results.pca_proj_matrix = pca_proj_matrix;

act_svm_results.TRloss


%% save active weight params
savefilename = 'wtparams_active.mat';
wtprms_active = [];

wtprms_active.tslen = tslen;
wtprms_active.svm_prms = act_svm_results;



save(savefilename, 'wtprms_active');


