%%  load data
load newdata;

dataset = 12;

trace_means = alldata(dataset).trace_mean;

%%  preprocess data

trace_means_in = trace_means';
[T,nroi] = size(trace_means_in);

norm_means = nan(size(trace_means_in));
norm_means_smooth = nan(size(trace_means_in));

x = 700;
p = 5;
e = 30;
b = 100;

filtersd = 4;
for t=1:T
    tic
    [norm_means(t,:), norm_means_smooth(t,:)] = preprocess_trace('F',trace_means_in, 'dFF', norm_means, 'timeIndex', t, ...
                                                                 'windowsize', x, 'percentile', p, 'epsilon', e, 'background', b, 'filtersd',filtersd);
                                                             toc
end


%%  visualization

trace_means_in;
norm_means;
norm_means_smooth;
nroi;


figure; 
hold on;
for roi=1:nroi
    plot(trace_means_in(:,roi), '-');
   % plot(norm_means(:,roi), '-');
   % pause;
end
title(sprintf('%s -- raw',alldata(dataset).name),'interpreter','none');

figure; 
hold on;
for roi=1:nroi
   % plot(trace_means_in(:,roi), '-', 'linewidth', 2);
    plot(norm_means(:,roi), '-');
   % pause;
end
title(sprintf('%s -- normalized',alldata(dataset).name),'interpreter','none');

figure; 
hold on;
for roi=1:nroi
   % plot(trace_means_in(:,roi), '-', 'linewidth', 2);
    plot(norm_means_smooth(:,roi), '-');
   % pause;
end
title(sprintf('%s -- smooth',alldata(dataset).name),'interpreter','none');

