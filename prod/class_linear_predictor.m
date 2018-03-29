%  
%  calls    fitclinear at its core
%
%
%  common lin_param = 
%        'Learner'  --  'svm' (default) | 'logistic'
%        'Regularization' --  'lasso' | 'ridge' (default)
%        'Cost'  --  ones(K) - eye(K)  (misclassification cost)  (default)
%        'Weights'  --  how much does each observation weight
%        'Lambda' -- (default = 1/n)
%        'OptimizeHyperparameters'  -- 'auto' (lambda and learner) | 'all' | {'lambda'}
%  for full documentation, see:
%     https://www.mathworks.com/help/stats/fitclinear.html
%
% **NOTE** IF PROVIDING 'Weights', 'Cost' or 'Prior', make sure that  'trian_balance' (1) == false
%
%  EXAMPLES:
%




classdef class_linear_predictor
    methods(Static)
  

        function lin_params_aug = augment_parameters_balance(lin_params)
            if isempty(lin_params) || ~strcmp(lin_params{1},'train_test_balance')
                lin_params_aug = {'train_test_balance',[false false], lin_params{:}};  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else
                lin_params_aug = lin_params;
            end
        end
        
        
        function [min_avgTSloss, all_TSloss, min_avgTRloss, all_TRloss, cv_ndx] = linpredict_cv(odat, labels, lin_params, cvndx)    
            if nargin<4 || isempty(cvndx)
                nfolds = 5;
                cvndx = cross_validate_split(lengthr(odat),nfolds);
            end
            lin_params = class_linear_predictor.augment_parameters_balance(lin_params);
            
            [min_avgTSloss, all_TSloss, min_avgTRloss, all_TRloss, cv_ndx] =  class_linear_predictor.do_cross_validate_linpredict(odat,labels, lin_params, cvndx);
        end
        
        function [wt_vec,lin_struct,TRloss,plabel] = linpredict(odat,labels, lin_params)
            lin_params = class_linear_predictor.augment_parameters_balance(lin_params);
            
            param_out = class_linear_predictor.train_fxn(odat,labels,[],lin_params); 
            lin_struct = param_out.linstruct;
            
            wt_vec = lin_struct.Beta;
% %             [pts,dim] = size(odat);
% %             a = lin_struct.Alpha;
% %             a_all = zeros(pts,1);
% %             a_all(lin_struct.IsSupportVector) = a;
% %             
% %             lbls = labels;
% %             lbls(lbls==0) = -1;
% %             
% %             
% %             wt_vec = sum(repmat(a_all,1,dim).*repmat(lbls,1,dim).*odat, 1);
% %             
            plabel = class_linear_predictor.pred_fxn(odat, param_out);
            
           to_test_balance = lin_params{2}(2);
           if to_test_balance
               loss_fxn_to_use = @class_linear_predictor.loss_fxn_balanced;
           else
               loss_fxn_to_use = @class_linear_predictor.loss_fxn;
           end
           
           
            TRloss = loss_fxn_to_use(plabel, labels)
        end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               
        
        
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function [min_avgloss, all_loss, min_avgTRloss, all_TRloss, cross_valid_ndx] = do_cross_validate_linpredict(x,y, prm, cross_valid_ndx)
           if nargin<4
               %cross_valid_ndx = 1:lengthr(x);
               cross_valid_ndx = cross_validate_split(lengthr(x), 5);
           end
           
           to_test_balance = prm{2}(2);
           if to_test_balance
               loss_fxn_to_use = @class_linear_predictor.loss_fxn_balanced;
           else
               loss_fxn_to_use = @class_linear_predictor.loss_fxn;
           end
           
           [min_param_setting, min_avgloss, all_loss, min_avgTRloss, all_TRloss] = ...
                 cross_validate(x, y, cross_valid_ndx, {[]}, prm, [], @class_linear_predictor.train_fxn, @class_linear_predictor.pred_fxn, loss_fxn_to_use);
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
               param_out.linstruct =  fitclinear(train_x,train_y,'Weights',wt_vec,train_params_in{:});
           else
               param_out.linstruct =  fitclinear(train_x,train_y,train_params_in{:});
           end
       end
       
       function pred_y = pred_fxn(pred_x, param_out, common_test_params)
           pred_y = predict(param_out.linstruct, pred_x);
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
