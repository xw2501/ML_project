classdef class_process_data
    methods(Static)
        
        function lblid = resolvelbl(lbl)
            lblid = [];
            if isempty(lbl); return; end;
            switch(lbl)
                case 'fwd'
                    lblid = 1;
                case 'weak_fwd'
                    lblid = 2;
                    
                case 'bkd'
                    lblid = 3;
                case 'weak_bkd'
                    lblid = 4;
                    
                case 'activity'
                    lblid = 5;
                case 'no_activity'
                    lblid = 6;
                    
                case 'no_fwdbkd'
                    lblid = 0;
                    
                otherwise
                    error('unrecognized lbl: %s',lbl);
            end
        end
        function lbldesc = get_lbldescr
            lbldesc = {'fwd' 'weak_fwd' 'bkd' 'weak_bkd' 'activity' 'no_activity'};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function alldata = load_and_prep_data(base_dir, data_dir, numrois)
            %%numrois = 6;
            %numrois = 18;
            %
            %base_dir = '/groups/branson/home/verman/data/work/projs/chen';
            %data_dir = sprintf('%s/data/%0dROIs',base_dir,numrois);
            
            
            dirs = dir(data_dir);
            
            alldata = [];
            
            
            for di = 1:length(dirs)
                dr = dirs(di);
                if ~dr.isdir; continue; end;
                if strcmp(dr.name,'.') || strcmp(dr.name,'..'); continue; end;
                
                ad = [];
                orig_data = [];
                
                dataloc = sprintf('%s/%s',data_dir,dr.name);
                try
                    ld_trace_dev = load(sprintf('%s/TraceDeviations.mat',dataloc));
                catch
                    warning(sprintf('File not found: %s/TraceDeviations.mat',dataloc));
                    ld_trace_dev = [];
                    ld_trace_dev.traceDeviations = [];
                end
                ld_trace_mean = load(sprintf('%s/TraceAverages.mat',dataloc));
                ld_anno = load(sprintf('%s/TraceAnnotations.mat',dataloc));
                try
                    ld_isactivity_anno = load(sprintf('%s/IsActivityAnnotations.mat',dataloc));
                catch
                    warning('IsActivityAnnotations.mat not found');
                    ld_isactivity_anno = [];
                    ld_isactivity_anno.IsActivityAnnotations = [];
                end
                
                orig_data.trace_dev = ld_trace_dev.traceDeviations;
                orig_data.trace_mean = ld_trace_mean.traceAverages;
                ad.name = dr.name;
                ad.orig_data = orig_data;
                
                if size(orig_data.trace_mean,1)~=numrois
                    ad.trace_mean = reshape(orig_data.trace_mean,[numrois size(orig_data.trace_mean,3)]);
                else
                    ad.trace_mean = orig_data.trace_mean;
                end
                ad.anno = [ld_anno.allAnnotations; ld_isactivity_anno.IsActivityAnnotations;];
                ad.anno_descr = class_process_data.get_lbldescr();
                alldata  = [alldata ad];
            end
            %alldata = class_process_data.normalize_data(alldata);
        end
        
        
        function alldata = normalize_data(alldata)
            
            stdev = 10; %frames
            filter = normpdf([-stdev*3:0]',0,stdev);
            filter = filter./sum(filter);
            
            for triali=1:length(alldata)
                ad = alldata(triali);
                dsz=size(ad.trace_mean);
                trace_normalized = nan(dsz);
                trace_mean = ad.trace_mean;
                for i=1:dsz(1)%*dsz(2)
                    %[a1,a2]=ind2sub(dsz(1),i);
                    t = trace_mean(i,:)';
                    b=nan(size(t));
                    n=nan(size(t));
                    for j=1:length(t)
                        %b(j) = min(t(1:j));
                        b(j) = prctile(t(1:j),1);
                        %%% filtered
                        tmp = t(1:j) - b(1:j);
                        curfilt = ffffilter(fixrange(end-j+1):end); curfilt = curfilt/sum(curfilt);
                        n(j) = sum(curfilt.*tmp(fixrange(-stdev*3+j):j));
                        % %                             %%%  unfiltered
                        % %                         n(j) = t(j)-b(j);
                    end
                    %trace_normalized(a1,a2,:) = n;
                    trace_normalized(i,:) = n;
                end
                alldata(triali).trace_normalized = trace_normalized;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [pos_ndxs,neg_ndxs,ign_ndxs] = get_ndxs(alldata,plbls,pin,pig, nlbls,nin,nig, ilbls,iin, display_result)
            % extracts lbl data s.t. plbls are labeled pos
            %                        nlbls are labeled neg
            %                        ilbls are ignored from BOTH pos and neg
            %
            %               pin/pig  tells exactly which offsets from start/end get included  (pin)
            %                                    and which offsets from start/end get ignored (pig)
            %
            %               same for nin/nig and iin/iig
            %                  (further details below)
            %
            % if pos is empty, complement of neg labels (ignoring the ignore range) will be considered positive
            % if neg is empty, complement of pos labels (ignoring the ignore range) will be considered negative
            %
            %   {p,n,i}{in} <==> contains two elements, which tells what is considered in range
            %   example  pin = {'st0' 'st5'} <==> pos_in is start+0 to start+5
            %            pin = {'st0' 'en0'} <==> pos_in is start+0 to end+0
            %            pin = {'st-5' 'en+5'} <==> pos_in is start-5 to end+5
            %
            %
            %
            %  example:  start of fwd wave:
            %               if you want pos label of first 3 timesteps of fwd wave
            %                           ignore a padding of 5 timesteps of around fwd wave (dont use them for neg)
            %                           ignore weak_fwd wave (dont use for neg)
            %                           rest is negative examples
            %
            %   get_ndxs(alldata,{'fwd'},  % only fwd is considered pos
            %                    {'st0 'st2'}, % include first three timesteps
            %                    {'st-5','en+5'}, % ignore a range of 5 timesteps around
            %                    {}, % neg is a complement of pos (no special neg lbls)
            %                    {},
            %                    {},
            %                    {'weak_fwd'}  % ignore weak_fwd
            %                    {'st0' 'en0'} % ignore the entire weak_fwd
            %           );
            %
            
            if nargin<10; display_result=false; end;
            
            
            pos_lblndxs = nan(length(plbls),1);
            for i = 1:length(plbls)
                pos_lblndxs(i) = class_process_data.resolvelbl(plbls{i});
            end
            
            neg_lblndxs = nan(length(nlbls),1);
            for i = 1:length(nlbls)
                neg_lblndxs(i) = class_process_data.resolvelbl(nlbls{i});
            end
            
            ign_lblndxs = nan(length(ilbls),1);
            for i = 1:length(ilbls)
                ign_lblndxs(i) = class_process_data.resolvelbl(ilbls{i});
            end
            
            assert(~(isempty(pos_lblndxs)&&isempty(neg_lblndxs)));  % both pos and neg cannot be empty
            assert(isempty(intersect(pos_lblndxs,neg_lblndxs))); % pos and neg should have nothing in common
            
            ntrials = length(alldata);
            pos_ndxs = cell(ntrials,1);
            neg_ndxs = cell(ntrials,1);
            ign_ndxs = cell(ntrials,1);
            
            for triali = 1:ntrials
                td = alldata(triali);
                rawy = td.anno;
                
                maxlen = size(td.trace_mean,2);
                posi = [];
                negi = [];
                igni = [];
                
                if isempty(pos_lblndxs); posi = 1:maxlen; end;
                if isempty(neg_lblndxs); negi = 1:maxlen; end;
                
                for yi = 1:size(rawy,1)
                    yst = fixrange(rawy(yi,1),maxlen);
                    yen = fixrange(rawy(yi,2),maxlen);
                    yid = rawy(yi,3);
                    
                    if ismember(yid, pos_lblndxs) % positive
                        [st,en] = class_process_data.resolvendxrange(yst,yen,pin);
                        posi = union(posi,st:en);                 % included as positive
                        negi = setdiff(negi, st:en);              % removed from negaitve
                        [st,en] = class_process_data.resolvendxrange(yst,yen,pig);
                        negi = setdiff(negi, st:en);              % posignores(padding around pos) are removed from negative
                    end
                    if ismember(yid, neg_lblndxs) % negative
                        [st,en] = class_process_data.resolvendxrange(yst,yen,nin);
                        negi = union(negi,st:en);                 % included as negative
                        posi = setdiff(posi, st:en);              % removed from positive
                        [st,en] = class_process_data.resolvendxrange(yst,yen,nig);
                        posi = setdiff(posi, st:en);              % negignores(padding around neg) are removed from positive
                    end
                    if ismember(yid, ign_lblndxs) % ignores
                        [st,en] = class_process_data.resolvendxrange(yst,yen,iin);
                        igni = union(igni,st:en);
                    end
                end
                posi = setdiff(posi, igni);              % ignores(padding around ing) are removed from
                negi = setdiff(negi, igni);              % BOTH positive and negative
                
                assert(isempty(intersect(posi,negi)));
                
                pos_ndxs{triali} = posi(:)';
                neg_ndxs{triali} = negi(:)';
                ign_ndxs{triali} = igni(:)';
            end
            
            if display_result
                for triali = 1:ntrials
                    td = alldata(triali);
                    rawy = td.anno;
                    maxlen = size(td.trace_mean,2);
                    
                    posi = pos_ndxs{triali};
                    negi = neg_ndxs{triali};
                    igni = ign_ndxs{triali};
                    
                    figure; hold on;
                    plot(posi,zeros(size(posi)), 'go');
                    plot(negi,zeros(size(negi)), 'ro');
                    plot(igni,zeros(size(igni)), 'bo');
                    
                    for yi = 1:size(rawy,1)
                        yst = fixrange(rawy(yi,1),maxlen);
                        yen = fixrange(rawy(yi,2),maxlen);
                        yid = rawy(yi,3);
                        
                        if ismember(yid, pos_lblndxs) % positive
                            plot([yst yen],[0 0], 'g-','linewidth',2);
                        end
                        if ismember(yid, neg_lblndxs) % negative
                            plot([yst yen],[0 0], 'r-','linewidth',2);
                        end
                        if ismember(yid, ign_lblndxs) % ignores
                            plot([yst yen],[0 0], 'b-','linewidth',2);
                        end
                    end
                    title(td.name,'interpreter','none');
                end
            end % display
            
        end % end function
        
        function [st,en] = resolvendxrange(yst,yen,yin)
            switch(yin{1}(1:2))
                case 'st'
                    base1 = yst;
                case 'en'
                    base1 = yen;
                otherwise
                    error('unresolved yin');
            end
            off1 = str2double(yin{1}(3:end));
            switch(yin{2}(1:2))
                case 'st'
                    base2 = yst;
                case 'en'
                    base2 = yen;
                otherwise
                    error('unresolved yin');
            end
            off2 = str2double(yin{2}(3:end));
            st = base1+off1;
            en=base2+off2;
        end
        
        
        function [x,y,trial_id,frame_id, traialnames, f_desc] = extract_training_generic(alldata,pos_ndxs,neg_ndxs,varargin)
            [tslen, datasource, point_density_pos, point_density_neg, min_intensity_thresh] = myparse(varargin, 'tslen',100, 'datasource', 'trace_normalized', 'point_density_pos', 1, 'point_density_neg', 1, 'min_intensity_thresh', -inf);
            x = [];
            y = [];
            trial_id = [];
            frame_id = [];
            
            trialnames = {};
            
            f_desc = [];
            
            for triali = 1:length(alldata)
                fprintf('processing trial %d/%d\n',triali,length(alldata));
                td = alldata(triali);
                traialnames = [trialnames; td.name];
                
                rawx = td.(datasource);
                %rawy = td.anno;
                
                low_intensity_ndx = find(max(rawx,[],1)<min_intensity_thresh);
                pos_ndxsi = pos_ndxs{triali};
                neg_ndxsi = neg_ndxs{triali};
                pos_ndxsi = pos_ndxsi(pos_ndxsi>tslen);
                neg_ndxsi = neg_ndxsi(neg_ndxsi>tslen);
                
                % remove min intensity points from training
                orig_sz = length(pos_ndxsi);
                pos_ndxsi = setdiff(pos_ndxsi,low_intensity_ndx);
                new_sz = length(pos_ndxsi);
                fprintf('\t fraction pos dropped from thresholding: %1.2f\n',(orig_sz - new_sz)/orig_sz);
                neg_ndxsi = setdiff(neg_ndxsi,low_intensity_ndx);
                
                %subsample pos
                pos_ndxsi = pos_ndxsi(:)';
                incl_idx_pos = rand(size(pos_ndxsi)) <= point_density_pos;
                pos_ndx_selected = pos_ndxsi(incl_idx_pos);
                
                %subsample neg
                neg_ndxsi = neg_ndxsi(:)';
                incl_idx_neg = rand(size(neg_ndxsi)) <= point_density_neg;
                neg_ndx_selected = neg_ndxsi(incl_idx_neg);
                
                
                trial_ndxs = [pos_ndx_selected neg_ndx_selected];
                trial_lbls = [ones(length(pos_ndx_selected),1); zeros(length(neg_ndx_selected),1)];
                
                % %                 trial_ndxs_all = [pos_ndxsi(:)' neg_ndxsi(:)'];
                % %                 trial_lbls_all = [ones(length(pos_ndxsi),1); zeros(length(neg_ndxsi),1)];
                % %
                % %                     % subsample
                % %                 incl_idx = rand(size(trial_ndxs_all)) <= point_density;
                % %                 trial_ndxs = trial_ndxs_all(incl_idx);
                % %                 trial_lbls = trial_lbls_all(incl_idx);
                % %
                
                x_triali = nan(length(trial_ndxs), numel(rawx(:,1:tslen)));
                
                i = 0;
                for cndx = trial_ndxs
                    i = i+1;
                    tmp = rawx(:,(cndx-tslen+1):cndx);
                    %x = [x; tmp(:)'];
                    x_triali(i,:) = tmp(:)';
                    
                    if isempty(f_desc)
                        assert(size(tmp,2)==tslen);
                        f_desc = arrayfun(@(i)(class_process_data.get_fdesc(i,size(tmp),tslen)), 1:numel(tmp),'uniformoutput',false);
                    end
                end
                x = [x; x_triali];
                y = [y; trial_lbls(:)];
                trial_id = [trial_id; repmat(triali,length(trial_lbls),1)];
                frame_id = [frame_id; trial_ndxs(:)];
            end
        end
        
        
        
        function [x,y,trial_id,frame_id, traialnames, f_desc] = extract_training_data_2lbl(alldata,poslbl,neglbl,varargin)
            [tslen, datasource, stride] = myparse(varargin, 'tslen',100, 'datasource', 'trace_mean','stride', 1);
            
            x = [];
            y = [];
            trial_id = [];
            frame_id = [];
            
            trialnames = {};
            
            f_desc = [];
            
            for triali = 1:length(alldata)
                
                td = alldata(triali);
                traialnames = [trialnames; td.name];
                
                rawx = td.(datasource);
                rawy = td.anno;
                nroi = size(rawx,1);
                
                pos_ndxs = class_process_data.getndxs_simple(rawy, poslbl, size(rawx,2),tslen,stride);
                neg_ndxs = class_process_data.getndxs_simple(rawy, neglbl, size(rawx,2),tslen,stride);
                
                trial_ndxs = [pos_ndxs(:)' neg_ndxs(:)'];
                trial_lbls = [ones(length(pos_ndxs),1); zeros(length(neg_ndxs),1)];
                xblock = nan(length(trial_ndxs),nroi*tslen);
                xbi = 1;
                for cndx = trial_ndxs
                    tmp = rawx(:,(cndx-tslen+1):cndx);
                    %x = [x; tmp(:)'];
                    xblock(xbi,:) = tmp(:)';
                    xbi=xbi+1;
                    if isempty(f_desc)
                        
                        assert(size(tmp,2)==tslen);
                        f_desc = arrayfun(@(i)(class_process_data.get_fdesc(i,size(tmp),tslen)), 1:numel(tmp),'uniformoutput',false);
                    end
                end
                x = [x; xblock];
                y = [y; trial_lbls(:)];
                trial_id = [trial_id; repmat(triali,length(trial_lbls),1)];
                frame_id = [frame_id; trial_ndxs(:)];
            end
        end
        
        function ndxs = getndxs_simple(rawy, lbl, maxlen,tslen,stride)
            if nargin<5; stride = 1; end;
            
            ndxs = [];
            lblndx = class_process_data.resolvelbl(lbl);
            ylbl = rawy(rawy(:,3)==lblndx,[1 2]);
            for yi = 1:size(ylbl,1)
                ndxs = [ndxs fixrange(ylbl(yi,1),maxlen):stride:fixrange(ylbl(yi,2),maxlen)];
            end
            ndxs = ndxs(ndxs>=tslen & ndxs<=maxlen);
        end
        
        
        function f_desc = get_fdesc(i,sz, tslen)
            [roii,timei] = ind2sub(sz,i);
            f_desc = sprintf('ROI:%d, TIME:%d',roii,timei-tslen);
        end
        
        
        function [x,y,trial_id,frame_id, traialnames] = extract_training_data(alldata,poslbl, varargin)
            [tslen, posneggap, notincl_lbl, posstep, negstep, datasource] = myparse(varargin, ...
                'tslen',100, 'posneggap',5, 'notincl_lbl',[], 'posstep',1,'negstep',5, 'datasource', 'trace_mean');
            
            x = [];
            y = [];
            trial_id = [];
            frame_id = [];
            
            trialnames = {};
            
            f_desc = [];
            
            for triali = 1:length(alldata)
                
                td = alldata(triali);
                traialnames = [trialnames; td.name];
                
                rawx = td.(datasource);
                
                rawy = td.anno;
                
                posid = class_process_data.resolvelbl(poslbl);
                notincl_id = class_process_data.resolvelbl(notincl_lbl);
                
                [pos_ndxs, neg_ndxs] = class_process_data.getndxs(rawy, posid, notincl_id, posneggap, posstep, negstep, tslen, size(rawx,3));
                %[pos_ndxs, neg_ndxs] = class_process_data.getstartndxs(rawy, posid, notincl_id, posneggap, posstep, negstep, tslen, size(rawx,3));
                
                trial_ndxs = [pos_ndxs(:)' neg_ndxs(:)'];
                trial_lbls = [ones(length(pos_ndxs),1); zeros(length(neg_ndxs),1)];
                for cndx = trial_ndxs
                    tmp = rawx(:,:, (cndx-tslen+1):cndx);
                    x = [x; tmp(:)'];
                    
                    
                    if isempty(f_desc)
                        assert(size(tmp,2)==2);
                        assert(size(tmp,3)==tslen);
                        f_desc = arrayfun(@(i)(class_process_data.get_feature_label(i,size(tmp))), 1:numel(tmp),'uniformoutput',false);
                    end
                end
                y = [y; trial_lbls(:)];
                trial_id = [trial_id; repmat(triali,length(trial_lbls),1)];
                frame_id = [frame_id; trial_ndxs(:)];
            end
        end
        
        function fname = get_feature_label(i, sz)
            [ai,di,ti] = ind2sub(sz,i);
            switch(di)
                case 1
                    d = 'left';
                case 2
                    d = 'right';
                otherwise
                    error('unknown direction');
            end
            fname = sprintf('A%d_%s_T%d',ai,d,ti-sz(3));
        end
        function [pos_ndxs, neg_ndxs] = getndxs(rawy, posid, notincl_id, posneggap, posstep,negstep, tslen, maxlen)
            warning('entire pos selected');
            
            pos_ndxs = [];
            neg_ndxs = [];
            if isempty(notincl_id); notincl_id = -999; end;
            
            allndxs = ones(1,maxlen);
            allndxs(1:tslen) = 0;
            
            posraw = rawy(rawy(:,end)==posid,:);
            if ~isempty(posraw)
                for i=1:size(posraw,1)
                    pr = posraw(i,:);
                    pos_ndxs = [pos_ndxs pr(1):posstep:pr(2)];
                    allndxs(  fixrange(pr(1)-posneggap,maxlen): fixrange(pr(2)+posneggap,maxlen) ) = 0;
                end
            end
            pos_ndxs = pos_ndxs(pos_ndxs>tslen);
            pos_ndxs = pos_ndxs(pos_ndxs<=maxlen);
            
            notinclraw = rawy(rawy(:,end)==notincl_id,:);
            if ~isempty(notinclraw)
                for i=1:size(notinclraw,1)
                    pr = notinclraw(i,:);
                    allndxs(  fixrange(pr(1)-posneggap,maxlen): fixrange(pr(2)+posneggap,maxlen) ) = 0;
                end
            end
            
            negraw = [find(diff([0 allndxs 0])==1)' find(diff([0 allndxs 0])==-1)'-1];
            if ~isempty(negraw)
                for i=1:size(negraw,1)
                    pr = negraw(i,:);
                    neg_ndxs = [neg_ndxs pr(1):negstep:pr(2)];
                end
            end
            
        end
        
        
        function [pos_ndxs, neg_ndxs] = getstartndxs(rawy, posid, notincl_id, posneggap, posstep,negstep, tslen, maxlen)
            warning('start ndxs only');
            pos_ndxs = [];
            neg_ndxs = [];
            
            posraw = rawy(rawy(:,end)==posid,:);
            if ~isempty(posraw)
                for i=1:size(posraw,1)
                    pr = posraw(i,:);
                    pos_ndxs = [pos_ndxs [0:5]+pr(1)];
                    neg_ndxs = [neg_ndxs [-5:-1 +6:+10]+pr(1)];
                end
            end
            pos_ndxs = pos_ndxs(pos_ndxs>tslen);
            pos_ndxs = pos_ndxs(pos_ndxs<=maxlen);
            
            neg_ndxs = neg_ndxs(neg_ndxs>tslen);
            neg_ndxs = neg_ndxs(neg_ndxs<=maxlen);
            
        end
        
        function [x,frame_id,f_desc] = all_instance_data(rawx,tslen)
            x = [];
            mx = size(rawx,2);
            frame_id = tslen:mx;
            f_desc = [];
            for cndx = frame_id
                tmp = rawx(:, (cndx-tslen+1):cndx);
                x = [x; tmp(:)'];
                if isempty(f_desc)
                    assert(size(tmp,2)==tslen);
                    f_desc = arrayfun(@(i)(class_process_data.get_fdesc(i,size(tmp),tslen)), 1:numel(tmp),'uniformoutput',false);
                    
                end
            end
        end
        
        
        function [train_X, train_Y, test_X,test_Y, trainidx] = train_test_split(all_X,all_Y, trial_id, frame_id,train_frac)
            trainidx = false(size(all_Y));
            ut = unique(trial_id);
            for ti = ut(:)'
                ti_idx = trial_id == ti;
                
                fi_idx = frame_id <= prctile(frame_id(ti_idx),train_frac*100);
                trainidx(ti_idx & fi_idx) = true;
            end
            
            train_X = all_X(trainidx,:);
            train_Y = all_Y(trainidx,:);
            test_X = all_X(~trainidx,:);
            test_Y = all_Y(~trainidx,:);
            
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end
end
