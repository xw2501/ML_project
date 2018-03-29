%  EXAMPLE svm parameter settings
%svm_cv_prm = {'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto'};
%svm_cv_prm = {'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto','Standardize',true};
%svm_cv_prm = {'KernelFunction', 'RBF'};
%svm_cv_prm = {'Standardize',true};
%svm_cv_prm = {};

% balanced svm ==>>  svm_cv_prm = {'train_test_balance', [true true], 'KernelFunction','polynomial','polynomialorder',1,'kernelscale','auto'};
%                         train_test_balance (if exists) must be the first parameter                     
classdef class_svm
    methods(Static)
        function new_lbls = remap_labels(old_lbls)
            % remaps arbitrary labels to 1 2 3 4 ...
            unq = unique(old_lbls);
            new_lbls = nan(size(old_lbls));
            for i=1:length(unq);
                u = unq(i);
                new_lbls(old_lbls==u) = i;
            end
        end
        
        function svm_params_aug = augment_parameters_balance(svm_params)
            if isempty(svm_params) || ~strcmp(svm_params{1},'train_test_balance')
                svm_params_aug = {'train_test_balance',[false false], svm_params{:}};  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                svm_params_aug = svm_params;
            end
        end
        
        
        function [min_avgTSloss, all_TSloss, min_avgTRloss, all_TRloss, cv_ndx] = svm_cv(odat, labels, svm_params, cvndx)    
            if nargin<4 || isempty(cvndx)
                nfolds = 5;
                cvndx = cross_validate_split(lengthr(odat),nfolds);
            end
            svm_params = class_svm.augment_parameters_balance(svm_params);
            
            [min_avgTSloss, all_TSloss, min_avgTRloss, all_TRloss, cv_ndx] =  class_svm.do_cross_validate_svm(odat,labels, svm_params, cvndx);
        end
        
        function [wt_vec,svm_struct,TRloss,plabel] = svm(odat,labels, svm_params)
            svm_params = class_svm.augment_parameters_balance(svm_params);
            
            param_out = class_svm.train_fxn(odat,labels,[],svm_params); 
            svm_struct = param_out.svmstruct;
            
            [pts,dim] = size(odat);
            a = svm_struct.Alpha;
            a_all = zeros(pts,1);
            a_all(svm_struct.IsSupportVector) = a;
            
            lbls = labels;
            lbls(lbls==0) = -1;
            
            
            wt_vec = sum(repmat(a_all,1,dim).*repmat(lbls,1,dim).*odat, 1);
            
            plabel = class_svm.pred_fxn(odat, param_out);
            
           to_test_balance = svm_params{2}(2);
           if to_test_balance
               loss_fxn_to_use = @class_svm.loss_fxn_balanced;
           else
               loss_fxn_to_use = @class_svm.loss_fxn;
           end
           
           
            TRloss = loss_fxn_to_use(plabel, labels)
        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       
% %        
% %        function [min_avgloss, all_loss, min_avgTRloss, all_TRloss, cross_valid_ndx] = do_cross_validate_svm(x,y, prm, cross_valid_ndx)
% %            [min_param_setting, min_avgloss, all_loss, min_avgTRloss, all_TRloss] = cross_validate(x, y, cross_valid_ndx, {[]}, prm, [], ...
% %                                                                                                   @class_lmnn.train_fxn, @class_lmnn.pred_fxn, @class_lmnn.loss_fxn);
% %        end
% %        
% %        function param_out = train_fxn(train_x, train_y, param_in, common_train_params)
% %            % common_train_params = {clf_k train_k Linit other_lmnn2_params}
% %            lmnn_input_params = common_train_params(2:end);
% %            
% %            [L] = lmnn2(train_x',train_y', lmnn_input_params{:});
% %            %L = my_itml(train_x,train_y,'initM', eye(size(train_x,2))*1e-6);
% %            
% %            param_out.clf_k = common_train_params{1};
% %            param_out.traink_k = common_train_params{2};
% %            
% %            param_out.L = L;
% %            param_out.train_x = train_x;
% %            param_out.train_y = train_y;
% %        end
% %        
% %        function pred_y = pred_fxn(pred_x, param_out, common_test_params)
% %            clf_k = param_out.clf_k;
% %            L = param_out.L;
% %            train_x = param_out.train_x;
% %            train_y = param_out.train_y;
% %            
% %            tfm_trainx = train_x * L';
% %            tfm_testx = pred_x * L';
% %            
% %            [idx,dst] = knnsearch(tfm_trainx, tfm_testx, 'K',clf_k+1);
% %            
% %            
% %            if sum(dst(:,1))==0; % training data
% %                subs = 2:clf_k+1;
% %            else % test data
% %             subs = 1:clf_k;
% %            end
% %            
% %            pred_y = mode(train_y(idx(:,subs)),2);
% %        end
% %        
% %        
% %        function loss = loss_fxn(pred_y, actual_y)
% %            loss = mean(abs(pred_y - actual_y));
% %        end

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function [min_avgloss, all_loss, min_avgTRloss, all_TRloss, cross_valid_ndx] = do_cross_validate_svm(x,y, prm, cross_valid_ndx)
           if nargin<4
               %cross_valid_ndx = 1:lengthr(x);
               cross_valid_ndx = cross_validate_split(lengthr(x), 5);
           end
           
           to_test_balance = prm{2}(2);
           if to_test_balance
               loss_fxn_to_use = @class_svm.loss_fxn_balanced;
           else
               loss_fxn_to_use = @class_svm.loss_fxn;
           end
           
           [min_param_setting, min_avgloss, all_loss, min_avgTRloss, all_TRloss] = ...
                 cross_validate(x, y, cross_valid_ndx, {[]}, prm, [], @class_svm.train_fxn, @class_svm.pred_fxn, loss_fxn_to_use);
       end
       
       function param_out = train_fxn(train_x, train_y, param_in, common_train_params)
           train_params_in = common_train_params(3:end);
           to_train_balance = common_train_params{2}(1);
           
           if to_train_balance
               tol = 1e-5;
                 % balancing code
               count_y = length(train_y);
               count_class0 = sum(train_y==0);
               count_class1 = sum(train_y==1);
               assert(count_class0+count_class1 == count_y);  % everything is either class 0 or class 1
               
               wt_vec = zeros(size(train_y));
               wt_vec(train_y==0) = count_y / count_class0 / 2;
               wt_vec(train_y==1) = count_y / count_class1 / 2;
               
               assert(count_class0 > 0 && count_class1>0);                            % at least one example from each class exists
               assert(abs(sum(wt_vec(train_y==0)) - sum(wt_vec(train_y==1))) < tol);  % weights are balanced
               assert(abs(sum(wt_vec) - count_y) < tol);                                        % weights sum up to number of points
               
                  % call balance weighted svm
               param_out.svmstruct =  fitcsvm(train_x,train_y,'Weights',wt_vec,train_params_in{:});
           else
               param_out.svmstruct =  fitcsvm(train_x,train_y,train_params_in{:});
           end
       end
       
       function pred_y = pred_fxn(pred_x, param_out, common_test_params)
           pred_y = predict(param_out.svmstruct, pred_x);
       end
       
       
       function loss = loss_fxn(pred_y, actual_y)
           loss = mean(abs(pred_y - actual_y));
       end

       function loss = loss_fxn_balanced(pred_y, actual_y)
           zero_one_err = abs(pred_y - actual_y);
           
              % balancing code
           tol = 1e-5;   
           count_y = length(actual_y);
           count_class0 = sum(actual_y==0);
           count_class1 = sum(actual_y==1);
           assert(count_class0+count_class1 == count_y);  % everything is either class 0 or class 1
           
           wt_vec = zeros(size(actual_y));
           wt_vec(actual_y==0) = 1 / 2 / count_class0; 
           wt_vec(actual_y==1) = 1 / 2 / count_class1;
           
           assert(count_class0 > 0 && count_class1>0);                                % at least one example from each class exists
           assert(abs( sum(wt_vec(actual_y==0)) - sum(wt_vec(actual_y==1)) ) < tol);  % equal weight is assigned
           assert(abs(sum(wt_vec) - 1) < tol);                                        % weights sum-up to one
           loss = sum(wt_vec .* zero_one_err);
       end


    end
end