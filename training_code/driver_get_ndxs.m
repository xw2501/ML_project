alldata;

%%

display_result = false;



plbls = {'fwd'};
pin = {'st0' 'en0'};
pig = {'st-5' 'en+5'};
nlbls = {}; 
nin = {}; 
nig = {};
ilbls = {'weak_fwd'};
iin = {'st0' 'en0'};

fwd_pos_ndxs = class_process_data.get_ndxs(alldata,plbls,pin,pig, nlbls,nin,nig, ilbls,iin, display_result);






plbls = {'bkd'};
pin = {'st0' 'en0'};
pig = {'st-5' 'en+5'};
nlbls = {}; 
nin = {}; 
nig = {};
ilbls = {'weak_bkd'};
iin = {'st0' 'en0'};

bkd_pos_ndxs = class_process_data.get_ndxs(alldata,plbls,pin,pig, nlbls,nin,nig, ilbls,iin, display_result);


%%
alldata;


for dsi = 1:length(alldata)
    ds = alldata(dsi);
    
    dff_smooth = ds.dFF_smooth';
    
    fw = fwd_pos_ndxs{dsi};
    bw = bkd_pos_ndxs{dsi};
    
    %fval = expand_feature_space(dff_smooth,[fw bw]);
    fval = expand_feature_space2(dff_smooth,[fw bw]);

    p = fval;
    x = length(fw);
    
    %p = diff(fval);
    %x = length(fw)-1;
    
    figure;
    hold on;
    plot(p,'.b');
    
    yl = ylim;
    xl = xlim;
    plot([x x],yl,'-r','linewidth',2);
    plot(xl,[5 5],'-k');
    plot(xl,[-5 -5],'-k');
    ylabel('feature value');
    
    title(sprintf('%s',ds.name),'interpreter','none');
 %   pause;
end

%%

