function [min_param_setting, min_avgloss, all_loss, min_avgTRloss, all_TRloss] = cross_validate(all_x, all_y, cross_valid_ndx, param_settings, common_train_params, common_test_params, train_fxn, pred_fxn, loss_fxn)
%
%  function [min_param_setting, min_avgloss, all_loss] = cross_validate(all_x, all_y, cross_valid_ndx, param_settings, common_params, train_fxn, pred_fxn, loss_fxn)
%
%  cross_validate over the various parameters settings specified in param_settings -- sweep over   1:length(param_settings)
%
%  min_param_setting is the param_settings{ndx}  which gives the lowest cross-validation error, wrt loss_fxn over splits induced by cross_valid_ndx
%
%
%  CALL STYLE
%
%  param_out = train_fxn(train_x, train_y, param_in, common_train_params);   %%%  param_in = one of the settings from param_settings,  ie param_settings{idx}
%  pred_y = pred_fxn(pred_x, param_out, common_test_params); 
%  loss = loss_fxn(pred_y, actual_y);                                               %%% must return a real number
%



num_folds = max(cross_valid_ndx);
num_param_settings = length(param_settings);


all_loss = nan(num_folds, num_param_settings);

all_loss_coll = nan(num_folds*num_param_settings,1);

all_TRloss = nan(num_folds, num_param_settings);
all_TRloss_coll = nan(num_folds*num_param_settings,1);

DISPLAY_PAR_PROGRESS = 0;
if DISPLAY_PAR_PROGRESS; parfor_progress(length(all_loss_coll)); end;
%
%  'parfor' loop
%
%parfor cidx = 1:num_folds*num_param_settings   % PARFOR   
for cidx = 1:num_folds*num_param_settings   % PARFOR
   
    [fold_num, param_num] = ind2sub([num_folds num_param_settings], cidx);
    
    valid_ndx = cross_valid_ndx==fold_num;
    trainx = all_x(~valid_ndx,:);
    trainy = all_y(~valid_ndx,:);
    validx = all_x(valid_ndx,:);
    validy = all_y(valid_ndx,:);
    
    param_in = param_settings{param_num};
    
    param_out = train_fxn(trainx, trainy, param_in, common_train_params);
    pred_y = pred_fxn(validx, param_out, common_test_params);
    
    all_loss_coll(cidx) = loss_fxn(pred_y, validy);

    pred_TRy = pred_fxn(trainx, param_out, common_test_params);
    all_TRloss_coll(cidx) = loss_fxn(pred_TRy, trainy);
    
    if DISPLAY_PAR_PROGRESS; parfor_progress; end;
end
if DISPLAY_PAR_PROGRESS; parfor_progress(0); end;



for cidx = 1:num_folds*num_param_settings    
    [fold_num, param_num] = ind2sub([num_folds num_param_settings], cidx);
    all_loss(fold_num,param_num) = all_loss_coll(cidx);
    all_TRloss(fold_num,param_num) = all_TRloss_coll(cidx);
end


avg_loss = mean(all_loss,1);  % avg loss over folds

[min_avgloss, min_avgloss_ndx] = min(avg_loss);
min_param_setting = param_settings{min_avgloss_ndx};

avg_TRloss = mean(all_TRloss,1);  % avg loss over folds

[min_avgTRloss, min_avgTRloss_ndx] = min(avg_TRloss);

end