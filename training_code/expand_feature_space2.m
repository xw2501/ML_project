function dff_smooth_expand = expand_feature_space2(dff_smooth, ndxs)

    nroi = size(dff_smooth,2);
    
    assert(nroi==18);
    
    lt = 1:9;
    rt = 10:18;  % agreed upon convention (see runTCPserver.m)
    
    
    
    out = arrayfun(@(ndx)(get_mx(dff_smooth,ndx,lt,rt)),ndxs); %,'uniformoutput',false);
    dff_smooth_expand = out;
end


function out = get_mx(dff_smooth,ndx,lt,rt)

TRAIL = 20;
ndx_stt = ndx - TRAIL;
ndx_end = ndx;

if ndx_stt<1; out=0; return; end;

[mx,mxndx] = max(dff_smooth(ndx_stt:ndx_end,:),[],1);

out = sum(sign(diff(mxndx(lt)))) + sum(sign(diff(mxndx(rt))));

end
