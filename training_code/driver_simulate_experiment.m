%% ensure stuff exists
alldata;
expid = 14;

%% do simulation & visualization


do_visualization = 1;
datasource = 'dFF_smooth';

if do_visualization;
    figure;
end

data_in = alldata(expid).(datasource)';
annos = alldata(expid).anno;
anno_descr = alldata(expid).anno_descr;
[T,nroi] = size(data_in);

online_pred_all = nan(5,T);

tic;
online_prediction('init', 1e5);
fprintf('INIT: '); toc;
%pause;
for t=1:T
    tic;
    [online_pred_all(:,t)] = online_prediction('data',data_in ,'timeIndex',t);
    fprintf('%d: ',t); 
    toc;
    
    if do_visualization
        is_fwd = online_pred_all(1,1:t);
        is_bkd = online_pred_all(3,1:t);
        is_active = online_pred_all(5,1:t);
        
        hold off;
        plot(max(data_in(1:t,:),[],2),'-b');
        hold on;

        plot(is_active, 'og');
        plot(is_fwd, '+r');
        plot(is_bkd, 'xb');
        
        legend('dFF-smooth','active','fwd','bkd');
        
        for ai = 1:size(annos,1)
            a = annos(ai,:);
            if a(1)>t; continue; end;
            clr = [];
            switch anno_descr{a(3)}
                case 'fwd'
                    clr = 'r';
                case 'weak_fwd'
                    clr = 'r';
                case 'bkd'
                    clr = 'b';
                case 'weak_bkd'
                    clr = 'b';
                    
                case 'activity'
                case 'no_activity'
            end
            if isempty(clr); continue; end;
            p = patch([a(1) a(2) a(2) a(1)],[0 0 1 1],clr);
            set(p,'facealpha',0.1, 'edgecolor','none');
            xl = xlim;
            xlim([2999 xl(2)]);
            
        end
        drawnow;
    end
end
