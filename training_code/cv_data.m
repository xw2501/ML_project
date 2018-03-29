function losses = cv_data(alldata)

trfxn = @tr_quick_fwd;
tsfxn = @ts_quick_fwd;

nfolds = length(alldata);

fld_range = 1:nfolds;
losses = cell(nfolds,1);
for foldi = fld_range
    ts_idx = fld_range==foldi;
    ts_alldata = alldata(ts_idx);
    tr_alldata = alldata(~ts_idx); 
    
    
    trprms = trfxn(tr_alldata);
    losses{foldi} = tsfxn(ts_alldata,trprms);
end

end


function trprms = tr_quick_fwd(tr_alldata)
    trprms = [];
end

function loss = ts_quick_fwd(ts_alldata, trprms)
    assert(length(ts_alldata)==1);
    ds = ts_alldata(1);
    dff_smooth = ds.dFF_smooth';
    T = size(dff_smooth,1);
    y =  expand_label(T, ds.anno, ds.anno_descr, 'fwd', 'weak_fwd');
    
    wtprms = [];
    wtprms.svm_prms.min_intensity_thresh = 0.015;
    predy = nan(1,T);
    fvals = nan(1,T);
    for t = 1:T
        fval = compute_special_feature(dff_smooth, t);
        fval_o1 = compute_special_feature(dff_smooth, t-1);
        fval_o2 = compute_special_feature(dff_smooth, t-2);
        fvals(t) = fval;
        
        predy(t) = classify_data(dff_smooth, t, wtprms, @quick_pred_fwd, {fval, fval_o1, fval_o2});
        
        %if y(t)==1
        %fprintf('%d -- %d\n', y(t), predy(t));
        %end
    end
    loss = compute_stats(y, predy);
f = figure; hold on; 
plot(y,'go'); 
plot(predy,'rx'); 
plot(fvals/5,'-b.');
title(ts_alldata.name,'interpreter','none'); legend('true-y','pred-y','fval');

%saveas(f,sprintf('fig_%s.fig',ts_alldata.name));
%close(f);
end

function loss = compute_stats(orig_y,orig_predy)
    ynanidx = isnan(orig_y);
    
        % y has no nans  (ie, no weak fwd/bkd anno);
    y = orig_y(~ynanidx);
    predy = orig_predy(~ynanidx);   
    
    tp = sum(predy==1 & y==1);  % true pos
    tn = sum(predy==0 & y==0);  % true neg
    fp = sum(predy==1 & y==0);  % false pos
    fn = sum(predy==0 & y==1);  % false neg
    rp = sum(isnan(predy) & y==1); % refrain prediction on pos
    rn = sum(isnan(predy) & y==0); % refrain prediction on neg
    
    pr = sum(predy==1 & isnan(y));  % pos prediction on refrained 
    nr = sum(predy==0 & isnan(y));  % neg prediction on refrained 
    rr = sum(isnan(predy) & isnan(y));  % refrain prediction on refrained
    
    loss = [tp, tn, fp, fn, rp, rn, pr, nr, rr];
end


function y = expand_label(len, anno, anno_descr, pos_anno, ig_anno)
    y = zeros(1,len);
    pndx = find(strcmp(anno_descr,pos_anno),1);
    indx = find(strcmp(anno_descr,ig_anno),1);
    
    a = anno(anno(:,3)==pndx,:);
    na = size(a,1);
    for i=1:na
        y(a(i,1):a(i,2)) = 1;
    end
    
    a = anno(anno(:,3)==indx,:);
    na = size(a,1);
    for i=1:na
        y(a(i,1):a(i,2)) = nan;
    end
end



function [result] = classify_data(dFF_smooth, nt, wtprms, quick_pred, inprms)
    result = NaN;

        % low dFF
    if max(abs(dFF_smooth(nt,:))) < wtprms.svm_prms.min_intensity_thresh;
        result = 0;
        return;
    end
  
        % clear signal
    tmpres = quick_pred(inprms);
    if ~isnan(tmpres); result = tmpres; return; end;
  
    if ~isfield(wtprms.svm_prms,'svm_struct'); return; end;
    
    tslen = wtprms.tslen;  
    if nt < tslen; return; end;    
    
        % unclear signal, checking further
    tmp = dFF_smooth((nt-tslen+1):nt, :)';  % remove baseline trend removed and filtered
    rawx = tmp(:)';
    
    normx = (rawx - wtprms.svm_prms.mn) ./ wtprms.svm_prms.st;
%    projx = normx * wtprms.pca_proj_matrix;

    result = class_svm.pred_fxn(normx, struct('svmstruct',wtprms.svm_prms.svm_struct)) * 0.5;
end





function is_fwd = quick_pred_fwd(inprm)
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
function is_bkw = quick_pred_bkd(fval)
    is_bkw = nan;
    if fval >= +5; is_bkw = 1; end
    if fval <= -5; is_bkw = 0; end
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
