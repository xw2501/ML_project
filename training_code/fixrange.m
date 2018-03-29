function [i] = fixrange(i,maxr,minr)
    if nargin<2; maxr=inf; end;
    if nargin<3; minr=1; end;
    if i>maxr; i=maxr; end;
    if i<minr; i=minr; end;
end