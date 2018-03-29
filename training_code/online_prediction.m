function [all_predictions] = online_prediction(varargin)
%
%   Example usage:
%
% 
%       >>   online_prediction_fwd_wave('init', buffersize);                                % call it once to initialize
%       >>   online_prediction_fwd_wave('data', roiBuffer, 'timeIndex', roiBufferTime);     %  call each time new data is collected
%
%
%
%   input:    'data'        roiBuffer_dFF_smooth      shape:  MAX_num_timepoints x (18)  [for 18 roi]
%             'timeIndex'   roiBufferTime             current time index
%
%   output:   5x1 vector, where each entry means:
%                   entry 1:    is_fwd_wave
%                   entry 2:    is_fwd_wave_start
%                   entry 3:    is_bkw_wave
%                   entry 4:    is_bkw_wave_start
%                   entry 5:    is_activity
%
%             each entry of the output vector can take three values:
%                            1  -- behavior is happening)
%                            0  -- behavior is NOT happening)  
%                           nan -- no prediction)
%
%
%



is_fwd_wave = nan; 
is_fwd_wave_start = nan;
is_bkw_wave = nan;
is_bkw_wave_start = nan;
is_active = nan;


all_predictions = [
                    is_fwd_wave;
                    is_fwd_wave_start;
                    is_bkw_wave;
                    is_bkw_wave_start;
                    is_active;    
                  ];
              

[to_init, dFF_smooth, nt] = myparse(varargin, 'init', 0, 'data', [], 'timeIndex', 0);

if to_init < 0 || floor(to_init)~=to_init; error('Initialization buffer size must be a positive integer'); end;
    
NROI = 18;

persistent wtparams_fwd;
persistent wtparams_bkd;
%persistent wtparams_act;
if to_init || isempty(wtparams_fwd); 
    wtparams_fwd = load('wtparams_fwd'); 
    wtparams_bkd = load('wtparams_bkw'); 
%    wtparams_act = load('wtparams_active');
    return;
end
wtprms_fwdwave = wtparams_fwd.wtprms_fwdwave;
wtprms_bkdwave = wtparams_bkd.wtprms_bkdwave;
%wtprms_active = wtparams_act.wtprms_active;


[nroi] = size(dFF_smooth,2);
if nroi~=NROI; error('input data must be of size:  MAX_buffer_size X %d ',NROI); end;


fval =  compute_special_feature(dFF_smooth, nt);
fval_o1 =  compute_special_feature(dFF_smooth, nt-1);
fval_o2 =  compute_special_feature(dFF_smooth, nt-2);


is_fwd_wave = classify_data(dFF_smooth, nt, wtprms_fwdwave, @quick_pred_fwd, {fval fval_o1 fval_o2});
is_bkw_wave = classify_data(dFF_smooth, nt, wtprms_bkdwave, @quick_pred_bkd, {fval fval_o1 fval_o2});
is_active   = quick_pred_active(dFF_smooth, nt); %classify_data(dFF_smooth, nt, wtprms_active, fval);


% CLEANUP

all_predictions = [
                    is_fwd_wave;
                    is_fwd_wave_start;
                    is_bkw_wave;
                    is_bkw_wave_start;
                    is_active;
                    ];


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result =  quick_pred_active(dFF_smooth, nt)
    result = 1;
    if max(abs(dFF_smooth(nt,:))) < 0.015
        result = 0.5;
    end
    if max(abs(dFF_smooth(nt,:))) < 0.01
        result = 0;
    end
end

function result = classify_data(dFF_smooth, nt, wtprms, quick_pred, inprms)
    result = NaN;

    tslen = wtprms.tslen;  

    if nt < tslen; return; end;

        % low dFF
    if max(abs(dFF_smooth(nt,:))) < wtprms.svm_prms.min_intensity_thresh;
        result = 0;
        return;
    end
  
        % clear signal
    tmpres = quick_pred(inprms);
    if ~isnan(tmpres); result = tmpres; return; end;
  
    if ~isfield(wtprms.svm_prms,'svm_struct'); return; end;
    
        % unclear signal, checking further
    tmp = dFF_smooth((nt-tslen+1):nt, :)';  % remove baseline trend removed and filtered
    rawx = tmp(:)';
    
    normx = (rawx - wtprms.svm_prms.mn) ./ wtprms.svm_prms.st;
%    projx = normx * wtprms.pca_proj_matrix;

    result = class_linear_predictor.pred_fxn(normx, struct('linstruct',wtprms.svm_prms.svm_struct)) * 0.5;
    %result = predict(normx, struct('svmstruct',wtprms.svm_prms.svm_struct)) * 0.5;
end





function is_fwd = quick_pred_fwd(inprm)
% %     is_fwd = nan;
% %     if fval <= -5; is_fwd = 1; end
% %     if fval >= +5; is_fwd = 0; end
% %     
fval = inprm{1};
fval_o1 = inprm{2};
fval_o2 = inprm{3};

    is_fwd = nan;
    if fval >= +5; is_fwd = 0; return; end
    if fval <= -7; 
        if fval_o1 < fval || fval_o2 < fval; return; end;
        if fval_o1 > fval && fval_o2 > fval_o1; is_fwd = 1; return ; end;    
    end
end
function is_bkw = quick_pred_bkd(inprm)
% %     is_bkw = nan;
% %     if fval >= +5; is_bkw = 1; end
% %     if fval <= -5; is_bkw = 0; end

fval = inprm{1};
fval_o1 = inprm{2};
fval_o2 = inprm{3};

    is_bkw = nan;
    if fval <= -5; is_bkw = 0; return; end
    if fval >= +7; 
        if fval_o1 > fval || fval_o2 > fval; return; end;
        if fval_o1 < fval && fval_o2 < fval_o1; is_bkw = 1; return ; end;    
    end

end


function fval = compute_special_feature(dff_smooth, nt)

    lt = 1:9;
    rt = 10:18;  % agreed upon convention (see runTCPserver.m)
    
    TRAIL = 20;
    ndx_stt = nt - TRAIL;
    ndx_end = nt;
    
    if ndx_stt<1; fval=0; return; end;
    
    [~,mxndx] = max(dff_smooth(ndx_stt:ndx_end,:),[],1);
    
    fval = sum(sign(diff(mxndx(lt)))) + sum(sign(diff(mxndx(rt))));
end

