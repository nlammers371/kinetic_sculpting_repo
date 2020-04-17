% Script to examine position-specific nucleus and activity flux trends over time 
close all
clear 
%%% set ID variables
project = 'eve7stripes_inf_2018_03_27_final';
DataPath = ['../../dat/' project];
FigPath = ['../../fig/experimental_system/' project '/stripe_dynamics/'];

mkdir(FigPath)
% load inference traces 
load([DataPath '\inference_traces_eve7stripes_inf_2018_02_20_dT20.mat']);
% load nuclei
load([DataPath '\inference_nuclei_eve7stripes_inf_2018_02_20_dT20.mat']);
%%% Define indexing vectors
% nc_stripe_vec = [nuclei_clean.stripe_id_vec_interp];
ts = [trace_struct_final.stripe_id_vec_interp];
ts = ts(~isnan(ts));
ts = round(3*ts)/3;
stripe_index = unique(round(ts(~isnan(ts))));
InterpGrid = round(trace_struct_final(1).InterpGrid);
%%% Generate Look-Up Tables for xp, yp and stripe ID
% traces
tr_x_array = NaN(length(InterpGrid),length(trace_struct_final));
tr_y_array = NaN(length(InterpGrid),length(trace_struct_final));
tr_s_array = NaN(length(InterpGrid),length(trace_struct_final));
tr_f_array = NaN(length(InterpGrid),length(trace_struct_final));
for i = 1:length(trace_struct_final)
    tt = round(trace_struct_final(i).time_interp);
    filter = ismember(InterpGrid,tt);
    ft = trace_struct_final(i).fluo_interp;
    xt = trace_struct_final(i).xPos_interp;
    yt = trace_struct_final(i).yPos_interp;
    st = trace_struct_final(i).stripe_id_vec_interp;
    tr_f_array(filter,i) = ft;
    tr_x_array(filter,i) = xt;
    tr_y_array(filter,i) = yt;
    tr_s_array(filter,i) = st;
end
tr_s_array = round(3*tr_s_array)/3;
% nuclei
nc_x_array = NaN(length(InterpGrid),length(nuclei_clean));
nc_y_array = NaN(length(InterpGrid),length(nuclei_clean));
nc_s_array = NaN(length(InterpGrid),length(nuclei_clean));

for i = 1:length(nuclei_clean)
    tn = round(nuclei_clean(i).time_interp);
    filter = ismember(InterpGrid,tn);    
    xn = nuclei_clean(i).xPos_interp;
    yn = nuclei_clean(i).yPos_interp;
    sn = nuclei_clean(i).stripe_id_vec_interp;
    
    nc_x_array(filter,i) = xn;
    nc_y_array(filter,i) = yn;
    nc_s_array(filter,i) = sn;
end
nc_s_array = round(3*nc_s_array)/3;
%%% --------------------- Fluo Flux Calculations -------------------------------- %%
n_steps = 15; % size of comparison window
inc = floor(n_steps/2);

% kernel for smoothing
kernel_radius = 6;
kernel_sigma = 3;
g_kernel = exp(-((-kernel_radius:kernel_radius)/(2*kernel_sigma)).^2);
g_kernel = g_kernel / sum(g_kernel);
% apply smoothing kernel to position and fluo arrays
ref_vec = ones(1,size(tr_x_array,1));
for i = 1:size(tr_f_array,2)
    % nuclei first
    nc_x = conv(nc_x_array(:,i)',g_kernel,'same')./conv(ref_vec,g_kernel,'same');
    nc_x_array(:,i) = nc_x;
    nc_y = conv(nc_y_array(:,i)',g_kernel,'same')./conv(ref_vec,g_kernel,'same');
    nc_y_array(:,i) = nc_y;
    % traces
    tr_x = conv(tr_x_array(:,i)',g_kernel,'same')./conv(ref_vec,g_kernel,'same');
    tr_x_array(:,i) = tr_x;
    tr_y = conv(tr_y_array(:,i)',g_kernel,'same')./conv(ref_vec,g_kernel,'same');
    tr_y_array(:,i) = tr_y;
    tr_f = conv(tr_f_array(:,i)',g_kernel,'same')./conv(ref_vec,g_kernel,'same');
    tr_f_array(:,i) = tr_f;
end
start_time = 30*60;
stop_time = 45*60;
calc_times = InterpGrid(InterpGrid>=start_time&InterpGrid<=stop_time);
calc_times = calc_times(1:end-inc);

tr_x_mean_flux_mat = NaN(length(calc_times),length(stripe_index));
tr_y_mean_flux_mat = NaN(length(calc_times),length(stripe_index));
t_ref_vec = InterpGrid(ceil(n_steps/2):end-floor(n_steps/2));

tr_x_diff = tr_x_array(n_steps:end,:) - tr_x_array(1:end-n_steps+1,:);
tr_y_diff = tr_y_array(n_steps:end,:) - tr_y_array(1:end-n_steps+1,:);
tr_f_wt = zeros(length(t_ref_vec),size(tr_f_array,2));

for i = inc+1:size(tr_f_array,1)-inc
    tr_f_wt(i,:) = sum(tr_f_array(i-inc:i+inc,:));
end
tr_s_ref = tr_s_array(ceil(n_steps/2):end-floor(n_steps/2),:);

for s = 1:length(stripe_index)
    stripe_id = stripe_index(s);
    for t = 1:length(calc_times)
        stripe_row = tr_s_ref(t_ref_vec==calc_times(t),:);
        x_start_vec = tr_x_array(InterpGrid==round(calc_times(t)-inc*20),stripe_row==stripe_id);
        x_stop_vec = tr_x_array(InterpGrid==round(calc_times(t)+inc*20),stripe_row==stripe_id);
        y_start_vec = tr_y_array(InterpGrid==round(calc_times(t)-inc*20),stripe_row==stripe_id);
        y_stop_vec = tr_y_array(InterpGrid==round(calc_times(t)+inc*20),stripe_row==stripe_id);
        f_start_vec = tr_f_array(InterpGrid==round(calc_times(t)-inc*20),stripe_row==stripe_id);
        f_stop_vec = tr_f_array(InterpGrid==round(calc_times(t)+inc*20),stripe_row==stripe_id);
        
        tr_x_mean_flux_mat(t,s) = nansum(x_stop_vec.*f_stop_vec)/nansum(f_stop_vec) - ...
                            nansum(x_start_vec.*f_start_vec)/nansum(f_start_vec);
        tr_y_mean_flux_mat(t,s) = nansum(y_stop_vec.*f_stop_vec)/nansum(f_stop_vec) - ...
                            nansum(y_start_vec.*f_start_vec)/nansum(f_start_vec);
    end
end

%%% ------------------- Nuclei Flux Calculations -------------------------------- %%
nc_x_flux_cell = cell(length(calc_times),length(stripe_index));
nc_x_mean_flux_mat = NaN(length(calc_times),length(stripe_index));
nc_y_flux_cell = cell(length(calc_times),length(stripe_index));
nc_y_mean_flux_mat = NaN(length(calc_times),length(stripe_index));

nc_x_diff = nc_x_array(n_steps:end,:) - nc_x_array(1:end-n_steps+1,:);
nc_y_diff = nc_y_array(n_steps:end,:) - nc_y_array(1:end-n_steps+1,:);
nc_s_ref = nc_s_array(ceil(n_steps/2):end-floor(n_steps/2),:);

for s = 1:length(stripe_index)
    stripe_id = stripe_index(s);
    for t = 1:length(calc_times)
        stripe_row = nc_s_ref(t_ref_vec==calc_times(t),:);
        x_deltas = nc_x_diff(t_ref_vec==calc_times(t),stripe_row==stripe_id);
        nc_x_flux_cell{t,s} = x_deltas(~isnan(x_deltas));
        nc_x_mean_flux_mat(t,s) = nanmean(x_deltas);
        y_deltas = nc_y_diff(t_ref_vec==calc_times(t),stripe_row==stripe_id);
        nc_y_flux_cell{t,s} = y_deltas(~isnan(y_deltas));
        nc_y_mean_flux_mat(t,s) = nanmean(y_deltas);
    end
end


%%% Fluo and Nuclei Flux Maps
cm = winter(length(calc_times));
cm2 = hot(length(calc_times))/1.2;

% inc = floor(7*128/7);
nf = conv(ones(1,size(nc_x_mean_flux_mat,1)),g_kernel,'same');
factor = n_steps/3;
nc_lim = 1.5;
nc_inc = (2*nc_lim)/10;
tr_lim_x = 15;
tr_inc_x = (2*tr_lim_x)/10;
tr_lim_y = 3;
tr_inc_y = (2*tr_lim_y)/10;
%%
for s = 1:7
    nc_flux_fig = figure;   
    hold on
    plot([-5 5],[0 0],'black')
    plot([0 0],[-5 5],'black')
    scatter(nc_x_mean_flux_mat(:,s)/factor,nc_y_mean_flux_mat(:,s)/factor,...
        linspace(5,40,length(calc_times)),cm,'filled','o');
    grid on        
    set(gca,'xtick',-nc_lim:nc_inc:nc_lim,'xticklabels',-nc_lim:nc_inc:nc_lim)
    set(gca,'ytick',-nc_lim:nc_inc:nc_lim,'yticklabels',-nc_lim:nc_inc:nc_lim)
    axis([-nc_lim nc_lim -nc_lim nc_lim])    
    xlabel('v_x (pixels per minute')
    ylabel('v_y (pixels per minute')
    title(['Mean Nuclear Drift: Stripe ' num2str(s) ' (' num2str(round(start_time/60))...
        ' - ' num2str(round(stop_time/60)) ')'])
    
    saveas(nc_flux_fig, [FigPath 'nc_dynamics_stripe' num2str(s) '.png'],'png')
    saveas(nc_flux_fig, [FigPath 'nc_dynamics_stripe' num2str(s) '.pdf'],'pdf')
    % particles
    tr_flux_fig = figure;
    hold on
    plot([-tr_lim_x tr_lim_x],[0 0],'black')
    plot([0 0],[-tr_lim_y tr_lim_y],'black')
    scatter(tr_x_mean_flux_mat(:,s)/factor,tr_y_mean_flux_mat(:,s)/factor,...
        linspace(5,40,length(calc_times)),cm2,...
        'filled','o');
    grid on        
    set(gca,'xtick',-tr_lim_x:tr_inc_x:tr_lim_x,'xticklabels',-tr_lim_x:tr_inc_x:tr_lim_x)
    set(gca,'ytick',-tr_lim_y:tr_inc_y:tr_lim_y,'yticklabels',-tr_lim_y:tr_inc_y:tr_lim_y)
%     axis([-nc_lim nc_lim -nc_lim nc_lim])    
    xlabel('v_x (pixels per minute)')
    ylabel('v_y (pixels per minute)')
    title(['Mean Drift in Transcriptional Activity: Stripe ' num2str(s) ' (' num2str(round(start_time/60))...
        ' - ' num2str(round(stop_time/60)) ')'])
    axis([-tr_lim_x tr_lim_x -tr_lim_y tr_lim_y])
    saveas(tr_flux_fig, [FigPath 'tr_dynamics_stripe' num2str(s) '.png'],'png')
    saveas(tr_flux_fig, [FigPath 'tr_dynamics_stripe' num2str(s) '.pdf'],'pdf')
    
    combined_flux_fig = figure;
    hold on
    plot([-tr_lim_x tr_lim_x],[0 0],'black')
    plot([0 0],[-tr_lim_y tr_lim_y],'black')
    scatter(tr_x_mean_flux_mat(:,s)/factor,tr_y_mean_flux_mat(:,s)/factor,...
        linspace(5,40,length(calc_times)),cm2,...
        'filled','o');
    scatter(nc_x_mean_flux_mat(:,s)/factor,nc_y_mean_flux_mat(:,s)/factor,...
        linspace(5,40,length(calc_times)),cm,'filled','o');
    grid on        
    set(gca,'xtick',-tr_lim_x:tr_inc_x:tr_lim_x,'xticklabels',-tr_lim_x:tr_inc_x:tr_lim_x)
    set(gca,'ytick',-tr_lim_y:tr_inc_y:tr_lim_y,'yticklabels',-tr_lim_y:tr_inc_y:tr_lim_y)        
    xlabel('v_x (pixels per minute)')
    ylabel('v_y (pixels per minute)')
    title(['Mean Drift Nuclei & Transcriptional Activity: Stripe ' num2str(s) ' (' num2str(round(start_time/60))...
        ' - ' num2str(round(stop_time/60)) ')'])
    axis([-tr_lim_x tr_lim_x -tr_lim_y tr_lim_y])
    saveas(combined_flux_fig, [FigPath 'combined_dynamics_stripe' num2str(s) '.png'],'png')
    saveas(combined_flux_fig, [FigPath 'combined_dynamics_stripe' num2str(s) '.pdf'],'pdf')
end        
        
%%
disp_fig = figure;
cm = jet(7);
hold on
plot([0 0],[1 7],'Color','Black')
scatter(mean(tr_x_mean_flux_mat)/factor*30,1:7,50,...
    cm,'filled')
plot([repelem(0,7) ; mean(tr_x_mean_flux_mat)/factor*30],[(1:7) ; (1:7)])

