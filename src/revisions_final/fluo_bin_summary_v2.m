% Script to generate summaries of inference results for fluo bin test
clear
close all
addpath('../utilities')
% Core parameters
K = 3; % State(s) to use for inference
w = 7; % Memory
t_start = 0; % start time for inference
dp_bootstrap = 1; % if 1 use bootstrap resampling at level of data points
project = 'revision_fluo_bins_v3';
ReadPath = '../../dat/revisions_final/';
FigPath = ['../../fig/revisions_final/' project '/'];
mkdir(FigPath);
WritePath = ['S:\Nick\Dropbox (Personal)\kinetic_scultping_inference_results\eve7stripes_data_v2\' project '/'];
mkdir(WritePath)
%%%% Stable Params (these rarely change) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
Tres = 20; % Time Resolution
min_dp = 10; % min length of traces to include

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([ReadPath 'inference_traces_w_stripe_ids.mat']);

alpha = trace_struct_final(1).alpha_frac*w; % Rise Time for MS2 Loops
read_dir = ['S:\Nick\Dropbox (Personal)\kinetic_scultping_inference_results\eve7stripes_data_v2\' ...
    project '\w' num2str(w) '_K' num2str(K) '_t' num2str(t_start) '_' project '\'];
% read_dir = ['../../out/revisions/' project '/w7_K3_revision_fluo_bins/'];
file_list = dir([read_dir '*.mat']);

inference_results = struct;
f_pass = 1;
for f = 1:numel(file_list)
    load([read_dir file_list(f).name]);
    fn = fieldnames(output);
    if numel(fn)>3
        for i = 1:numel(fn)
            inference_results(f_pass).(fn{i}) = output.(fn{i});
        end
        f_pass = f_pass + 1
    end
end

fluo_bin_vec = [inference_results.fluo_bin];
fluo_val_vec = [inference_results.fluo_val];
fluo_index = unique(fluo_bin_vec);
stripe_id_vec = [inference_results.Stripe];
stripe_index = unique(stripe_id_vec);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Checks to ensure validity of fluo bin identifiers %%%%%%%%%
tr_particle_index = [trace_struct_final.ParticleID];
% tr_stripe_index = [trace_struct_final.mike_stripe_id];
stripe_cx_flags = false(size(inference_results));
mean_fluo_vec = NaN(size(inference_results));
for i = 1:numel(inference_results)
    % inference info
    particle_ids = inference_results(i).particle_ids;
    stripe_id = stripe_id_vec(i);
%     fluo_id = fluo_bin_vec(i);
%     fluo_ids = tr_fluo_id(ismember(tr_particle_id,particle_ids));
    % cross-reference with trace characteristics
    stripe_cx_vec = NaN(size(particle_ids));
    fluo_cx_vec = NaN(size(particle_ids));
    for p = 1:length(particle_ids)
        p_ft = tr_particle_index == particle_ids(p);
        if t_start == 0
            time = trace_struct_final(p_ft).time_interp; 
            fluo = trace_struct_final(p_ft).fluo_interp; 
        else        
            time = trace_struct_final(p_ft).time_interp(w+1:end); 
            fluo = trace_struct_final(p_ft).fluo_interp(w+1:end);  
        end
        time_ft = time  >= t_start; 
        stripe_cx_vec(p) = mode(round(trace_struct_final(p_ft).mike_stripe_id(time_ft)));           
        fluo_cx_vec(p) = nanmean(fluo(time_ft));
    end	
    mean_fluo_vec(i) = nanmean(fluo_cx_vec);
    stripe_cx_flags(i) = all(stripe_cx_vec==stripe_id) || stripe_id == 0;
end

if ~all(stripe_cx_flags)
    error('inconsistency with stripe assignments')
end

close all
cm = brewermap(numel(stripe_index),'Set2');
for s = stripe_index
    fluo_cx_fig = figure;
    scatter(fluo_bin_vec(stripe_id_vec==s),mean_fluo_vec(stripe_id_vec==s),...
        'MarkerFaceColor',cm(s+1,:),'MarkerEdgeColor','black','MarkerEdgeAlpha',.2,'MarkerFaceAlpha',.3)
    xlabel('fluorescence bin')
    ylabel('inference group mean fluorescence (au)')
    set(gca,'Fontsize',14)
    xlim([0 10])
    ylim([0 2e5])
    grid on
    saveas(fluo_cx_fig,[FigPath 'data_integrity_check_stripe' num2str(s) '.png'])
end
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% iterate through results and generate summary stats %%%%%%%%%%%
hmm_results = struct;
transfer_vars = {'r','noise','soft_struct','Stripe','fluo_bin','fluo_val',...
    'particle_ids','N','w','alpha','deltaT'};
r_error = 0;
for i = 1:numel(inference_results)
    for j = 1:numel(transfer_vars)
        hmm_results(i).(transfer_vars{j}) = inference_results(i).(transfer_vars{j});
    end
    [r_sorted ,si] = sort(hmm_results(i).r);  
    % get transition rates 
    A = inference_results(i).A_mat(si,si);    
    dT = inference_results(i).deltaT;
    R = logm(A) / dT; 
    
%     % get occupancy
%     [V, D] = eigs(A);
%     [~, mi] = max(diag(D));
%     ss = V(:,mi);
%     ss = ss / sum(ss);                             
        
    % initialize fields
    hmm_results(i).r = NaN;
    hmm_results(i).A_mat = NaN;
    hmm_results(i).R = NaN;
    hmm_results(i).kon = NaN;
    hmm_results(i).kon_old = NaN;
    hmm_results(i).koff = NaN;
    hmm_results(i).koff = NaN;
    hmm_results(i).cycle_time = NaN;
    hmm_results(i).occupancy = NaN;
    hmm_results(i).r_eff = NaN;    
    hmm_results(i).r = NaN(K,1);
    if ~isreal(R) || sum(R(:)<0) > K
        out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
        R = out.R_out;     
    end
    % check for problematic init rate values
    r_error = r_error + 1*(abs(min(r_sorted)) > .05*max(r_sorted));
    
    % if inference results meet QC criteria, record
    if isreal(R) && sum(R(:)<0) == K && min(r_sorted) > -.05*max(r_sorted) 
        hmm_results(i).r = r_sorted;
        hmm_results(i).A_mat = A;
        hmm_results(i).R = R;
        
        % calculate effective 2 state parameters
        [V, D] = eigs(R);
        [~, mi] = max(diag(D));
        ss_vec = V(:,mi);
        ss_vec = ss_vec / sum(ss_vec); 
        
        hmm_results(i).kon = -R(1,1) ;
        hmm_results(i).kon_old = R(2,1) ;
        hmm_results(i).koff = -R(1,1)/(1/ss_vec(1) - 1);
        hmm_results(i).koff_old = R(1,2) * ss_vec(2) / sum(ss_vec(2:3));                
        hmm_results(i).r_eff = sum(r_sorted(2:3).*ss_vec(2:3))/sum(ss_vec(2:3));               
        hmm_results(i).occupancy = ss_vec;           
    end
end

% generate summary table to write to csv
save_prefix = ['w' num2str(w) '_K' num2str(K) '_t' num2str(t_start) '_'];
r_eff_vec = [hmm_results.r_eff];
kon_eff_vec = [hmm_results.kon];
koff_eff_vec = [hmm_results.koff];
r_mat = reshape([hmm_results.r]',numel(kon_eff_vec),K);

ft = ~isnan(kon_eff_vec);

inference_id_vec = 1:numel(hmm_results);
results_table = array2table([inference_id_vec' fluo_bin_vec' stripe_id_vec' kon_eff_vec' koff_eff_vec' ...
        r_mat r_eff_vec'], 'VariableNames',{'inf_id', 'fluo_bin', 'stripe_id',...
        'kon','koff','r1','r2','r3','r_eff'});
% save
writetable(results_table,[WritePath save_prefix 'fluo_bin_results_full_final.csv'])
save([WritePath '\' save_prefix 'fluo_bin_results_full_final.mat'],'hmm_results')

%%% Generate and save list of particle IDs corresponding to each inference
particle_id_cell = cell(1,numel(inference_results));
for i = 1:numel(inference_results)
    particle_id_cell{i} = inference_results(i).particle_ids;
end
save([WritePath save_prefix 'fluo_bin_particle_ids.mat'],'particle_id_cell')


% flag results with problematic r values for exclusion
r_err_flag = 0.05*abs(results_table.r1) > results_table.r3;

% make condensed summary
fluo_index = unique(fluo_bin_vec);
stripe_index = unique(stripe_id_vec);
fluo_val_index = unique([inference_results.fluo_val]);
iter = 1;
for f = 1:numel(fluo_index)
    for s = 1:numel(stripe_index)
        ft = fluo_bin_vec==fluo_index(f)&stripe_id_vec==stripe_index(s)&~isnan(kon_eff_vec)&~r_err_flag';   
        ft_raw = fluo_bin_vec==fluo_index(f)&stripe_id_vec==stripe_index(s);
        % obtain param vectors
        r_vec = r_eff_vec(ft);
        koff_vec = koff_eff_vec(ft);
        kon_vec = kon_eff_vec(ft);
        % flag outlier values
        [~,TF_r] = rmoutliers(r_vec);
        [~,TF_kon] = rmoutliers(kon_vec);
        [~,TF_koff] = rmoutliers(koff_vec);
        outlier_ft = ~(TF_r|TF_kon|TF_koff);
        % calculate statistics 
        r_med = 60*nanmean(r_vec(outlier_ft));
        r_ste = 60*nanstd(r_vec(outlier_ft));
        
        kon_med = 60*nanmean(kon_vec(outlier_ft));
        kon_ste = 60*nanstd(kon_vec(outlier_ft));
        
        koff_med = 60*nanmean(koff_vec(outlier_ft));
        koff_ste = 60*nanstd(koff_vec(outlier_ft));
        
        if iter == 1
            summary_mat = [fluo_index(f) fluo_val_index(f) stripe_index(s) sum(ft_raw) sum(ft) kon_med kon_ste koff_med koff_ste r_med r_ste];
        else
            summary_mat = vertcat(summary_mat,[fluo_index(f) fluo_val_index(f) stripe_index(s) sum(ft_raw) sum(ft) kon_med kon_ste koff_med koff_ste r_med r_ste]);
        end
        iter = iter + 1;
    end
end

summary_table = array2table(summary_mat(:,[1:3 6:11]), 'VariableNames',{'fluo_bin', 'mean_fluo', 'stripe_id',...
        'kon','kon_err','koff','koff_err','r','r_err'});
writetable(summary_table,[WritePath save_prefix 'fluo_hmm_results_final.csv'])
