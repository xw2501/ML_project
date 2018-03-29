function dff_smooth_expand = expand_feature_space(dff_smooth, ndxs)

    nroi = size(dff_smooth,2);
    
    assert(nroi==18);
    
   % left = [[1:2:nroi]
    
    lt = 1:9;
    rt = 10:18;  % agreed upon convention (see runTCPserver.m)
    
    %ndxs = 100:size(dff_smooth,1);
    
    
    for ti = 0:10
        lt_ftrs = arrayfun( @(roi)(make_ftr_space(roi,lt,dff_smooth, ndxs-ti)), 1:length(lt),'uniformoutput',false);
        rt_ftrs = arrayfun( @(roi)(make_ftr_space(roi,rt,dff_smooth, ndxs-ti)), 1:length(rt),'uniformoutput',false);

        coll = sum(cell2mat(lt_ftrs)+cell2mat(rt_ftrs),2);
        if ti==0
            f = coll;
        else
            f = f+coll;
        end
        
      
    end
    dff_smooth_expand = f;
    
end

function x_ftr = make_ftr_space(roii,dirroi,x, ndxs)
    nroi_trail_max = 5;
    nftr = min([nroi_trail_max roii-1]);
    if nftr<=0; x_ftr = []; return; end;
    
    x_ftr = nan(length(ndxs),nftr);
    
    for i=1:nftr
        x_ftr(:,i) = x(ndxs,dirroi(roii)) - x(ndxs,dirroi(roii-i));
%fprintf('idx:%d %d -- actual_roi: %d %d\n',roii, roii-i, dirroi(roii), dirroi(roii-i));        
    end
end

% % function x_ftr = make_ftr_time(time_trail_max,x,ndxs)
% % % assuming ndxs-ti will never be <= 0
% % 
% %     x_ftr = nan(length(ndxs),
% %     for ti = 1:time_trail_max
% %         
% %     end
% % end