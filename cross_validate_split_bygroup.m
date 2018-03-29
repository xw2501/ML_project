function ndx = cross_validate_split_bygroup(npts, num_folds, groupid)
    if isempty(npts)
        npts = length(groupid);
    end
    
    assert(npts == length(groupid));
    
    uniqgrps = unique(groupid);
    
    ngrps = length(uniqgrps);
    
    prm = randperm(ngrps);
    
    spl = round(linspace(1,ngrps, num_folds+1));
    
    ndx_grp = nan(ngrps,1);
    for i=1:num_folds
        if i==1
            ndx_grp(prm(spl(i):spl(i+1))) = i;
        else
            ndx_grp(prm(spl(i)+1:spl(i+1))) = i;
        end
    end

    
    ndx = nan(npts,1);
    
    for i=1:ngrps
        ndx(groupid==uniqgrps(i)) = ndx_grp(i);
    end
    
    if sum(isnan(ndx))>0
        error('unassigned split');
    end
end