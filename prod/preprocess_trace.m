function [dFt, dFts] = preprocess_trace(varargin)
%
%     data     MAX_T x NROI
%

[F,dF, t,x,p,e,b, filtersd] = myparse(varargin,'F', [], 'dFF', [], 'timeIndex',0, 'windowSize',0, 'percentile',0, 'background',0, 'epsilon',0, 'filterSD', 0);

nroi = size(F,2);

F0t = prctile(F(fixrange(t-x):t,:),p,1);

dFt = (F(t,:) - F0t) ./ (F0t - b + e);

%dFt(dFt<0) = 0;

dFts = dFt;

if filtersd<=0; return; end;

filter = getfilter(filtersd);
flen = length(filter);

for i=1:nroi
    tmp = [dF(fixrange(t-flen+1):t-1,i); dFt(i)];
    curfilt = filter(fixrange(end-t+1):end); curfilt = curfilt/sum(curfilt);
    dFts(i) = sum(curfilt.*tmp); 
end

end



function filter = getfilter(stdev)
    %stdev = 2; %frames
    filter = normpdf([-stdev*3:0]',0,stdev);
    filter = filter./sum(filter);
end

