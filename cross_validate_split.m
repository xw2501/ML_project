function ndx = cross_validate_split(npts, num_folds)
    prm = randperm(npts);
    
    spl = round(linspace(1,npts, num_folds+1));
    
    ndx = nan(npts,1);
    for i=1:num_folds
        if i==1
            ndx(prm(spl(i):spl(i+1))) = i;
        else
            ndx(prm(spl(i)+1:spl(i+1))) = i;
        end
    end
    if sum(isnan(ndx))>0
        error('unassigned split');
    end
end