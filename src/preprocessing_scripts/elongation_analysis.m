addpath('../utilities');
%------------------Define LoadPath Variables------------------------------%
datatype = 'weka';
project = 'eve7stripes_inf_2017_09_25';
Tres = 20;
%Path to raw data
datapath = ['../../dat/' project '/'];
start_time = 0;
stop_time = 60;

% generate save names
dataname_int = ['inference_traces_t' num2str(Tres) '_' project '.mat'];

% dataname_raw = 'raw_traces_s0_60_mHMMeve2_weka_inf_2017_09_25';

date_str = '5_partitions_test';
autopath = ['../../fig/experimental_system/' project '/preprocessing/ElongationTime/'];
if exist(autopath) ~= 7
    mkdir(autopath);
end
stripe_range = 1:7;
n_sets = length(stripe_range);
lags = 40;
% Load data 
load([datapath dataname_int]);


set_vec = unique([interp_struct.setID]);
%% By Region
a_avg_bin = NaN(lags+1,length(stripe_range));
a_std_bin = NaN(lags+1,length(stripe_range));
dd_avg_bin = NaN(lags-1,length(stripe_range));
dd_std_bin = NaN(lags-1,length(stripe_range));
for b = 1:length(stripe_range)
    s = stripe_range(b);
    bin_struct = interp_struct([interp_struct.stripe_id]==s);
    trace_array = zeros(20,length(bin_struct));
    for i = 1:length(bin_struct)
        trace = bin_struct(i).fluo;        
        trace_array(1:length(trace),i) = trace;
    end
    [am, as, ddm, dds] = weighted_autocorrelation(trace_array,lags,1);
    a_avg_bin(:,b) = am;
    a_std_bin(:,b) = as;
    dd_avg_bin(:,b) = ddm;
    dd_std_bin(:,b) = dds;
end

t_vec = 1:lags+1;
auto_fig = figure;
errorbar(repmat(t_vec',1,length(stripe_range)),a_avg_bin,a_std_bin,'LineWidth',1.5);
grid on
legend('Stripe 1', 'Stripe 2', 'Stripe 3', 'Stripe 4', 'Stripe 5', 'Stripe 6', 'Stripe 7')
title('Mean AutoCorrelation By Stripe');
saveas(auto_fig,[autopath 'mean_auto_by_bin_t' num2str(round(Tres)) '.png'],'png')

dd_fig = figure;
hold on
errorbar(repmat(t_vec(1:end-2)',1,length(stripe_range)),dd_avg_bin,dd_std_bin,'LineWidth',1.5);
grid on
legend('Stripe 1', 'Stripe 2', 'Stripe 3', 'Stripe 4', 'Stripe 5', 'Stripe 6', 'Stripe 7')
title('Mean Second Derivative of AutoCorrelation By Stripe');
saveas(dd_fig,[autopath 'mean_dd_by_bin_t' num2str(round(Tres)) '.png'],'png')
%% By Set 
legend_string = [];
lags = 40;
a_avg_set = NaN(lags+1,length(set_vec));
a_std_set = NaN(lags+1,length(set_vec));
dd_avg_set = NaN(lags-1,length(set_vec));
dd_std_set = NaN(lags-1,length(set_vec));
for b = 1:length(set_vec)
    set = set_vec(b);
    legend_string = [legend_string {['Set ' num2str(set)]}];
    set_struct = interp_struct([interp_struct.setID]==set);
    trace_array = zeros(20,length(set_struct));
    for i = 1:length(set_struct)
        trace = set_struct(i).fluo;                   
        trace_array(1:length(trace),i) = trace;
    end
    [am, as, ddm, dds] = weighted_autocorrelation(trace_array,lags,1);
    a_avg_set(:,b) = am;
    a_std_set(:,b) = as;
    dd_avg_set(:,b) = ddm;
    dd_std_set(:,b) = dds;
end

t_vec = 1:lags+1;
auto_fig = figure;
errorbar(repmat(t_vec',1,length(set_vec)),a_avg_set,a_std_set,'LineWidth',1);
grid on
legend(legend_string{:})
title('Mean AutoCorrelation By Set');
saveas(auto_fig,[autopath 'mean_auto_by_set_t' num2str(round(Tres)) '.png'],'png')

dd_fig = figure;
hold on
errorbar(repmat(t_vec(1:end-2)',1,length(set_vec)),dd_avg_set,dd_std_set,'LineWidth',1);
grid on
legend(legend_string{:})
title('Mean Second Derivative of AutoCorrelation By Set');
saveas(dd_fig,[autopath 'mean_dd_by_set_t' num2str(round(Tres)) '.png'],'png')


%% Full Set
lags = 40;
trace_array = zeros(20,length(interp_struct));
for i = 1:length(interp_struct)
    trace = interp_struct(i).fluo;    
    trace_array(1:length(trace),i) = trace;
end
[am, as, ddm, dds] = weighted_autocorrelation(trace_array,lags,1);
a_avg_set = am;
a_std_set = as;
dd_avg_set = ddm;
dd_std_set = dds;

t_vec = 1:lags+1;
auto_fig = figure;
errorbar(t_vec,a_avg_set,a_std_set,'LineWidth',1.5);
grid on
title('Mean AutoCorrelation All Sets');
saveas(auto_fig,[autopath 'mean_auto_all_set_t' num2str(round(Tres)) '.png'],'png')

dd_fig = figure;
hold on
errorbar(t_vec(1:end-2),dd_avg_set,dd_std_set,'LineWidth',1.5);
grid on
title('Mean Second Derivative of AutoCorrelation All Sets');
saveas(dd_fig,[autopath 'mean_dd_all_set_t' num2str(round(Tres)) '.png'],'png')

%% Time Trends in Mean Profile
stop_time = 60;

for i = 1:length(stripe_range)   
    
    bin_traces = interp_struct([interp_struct.stripe_id]==i);
    bin_traces = bin_traces([bin_traces.stripe_sub_id]==0);
    sets = unique([bin_traces.setID]);
    temporal_fig = figure('Position',[0 0 1024 512]);
    hold on
    time_vec = 1:stop_time;
    f_avg = zeros(1,length(time_vec));
    f_std = zeros(1,length(time_vec));
    legend_string = {};
    for j = 1:length(set_vec)+1
        if j <= length(set_vec)            
            set_struct = bin_traces([bin_traces.setID] == set_vec(j));
            lw = 1;            
        else            
            set_struct = bin_traces;
            lw = 2;
        end
        
        if isempty(set_struct)
            continue
        end
        time_list = ceil([set_struct.time]/60);
        fluo = [set_struct.fluo];
        for t = time_vec
            f_avg(t) = mean(fluo(ismember(time_list,[t-1:t+1])));
            f_std(t) = std(fluo(ismember(time_list,[t-1:t+1])));
        end
        if j <= length(set_vec)
            plot(time_vec,f_avg,'-o','LineWidth',1) 
            legend_string = {legend_string{:} ['Set ' num2str(set_vec(j))]};
        else
            plot(time_vec, f_avg, 'black', 'LineWidth',2);
            legend_string = {legend_string{:} 'Combined'};
        end
        
    end
    title(['Average Fluorescence Over Time (3 min moving avg), Stripe ' num2str(i)]);
    legend(legend_string{:},'Location','northwest')
    grid on    
    saveas(temporal_fig,[autopath 'time_trends_stripe' num2str(i) '.png'],'png')
end
close all