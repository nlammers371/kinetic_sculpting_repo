% Script to generate summaries of inference results for fluo bin test
clear
close all
addpath('../utilities')
% Core parameters
K = 3; % State(s) to use for inference
w = 7; % Memory
dp_bootstrap = 1; % if 1 use bootstrap resampling at level of data points
project = 'revision_fluo_bins';
ReadPath = '../../dat/revisions/';
FigPath = ['../../fig/revisions/' project '/'];
mkdir(FigPath);
WritePath = ['../../out/revisions/' project '/'];
mkdir(WritePath)
%%%% Stable Params (these rarely change) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
Tres = 20; % Time Resolution
min_dp = 10; % min length of traces to include

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([ReadPath 'fluo_inf_struct.mat']);
%%%
alpha = fluo_inf_struct(1).alpha_frac*w; % Rise Time for MS2 Loops

read_dir = ['../../out/revisions/' project '/w7_K3_revision_fluo_bins/'];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Checks to ensure validity of fluo bin identifiers %%%%%%%%%
tr_fluo_id = [fluo_inf_struct.FluoBin];
tr_particle_id = [fluo_inf_struct.ParticleID];
for i = 1:numel(inference_results)
    particle_ids = inference_results(i).particle_ids;
    fluo_ids = tr_fluo_id(ismember(tr_particle_id,particle_ids));
    fID = unique(fluo_ids);
    if numel(fID) > 1
        error('multiple fluo identifiers')
    else
        inference_results(i).fluo_bin = fID;
    end
end

%%% make indexing vectors
fluo_bin_vec = [inference_results.fluo_bin];
fluo_index = unique(fluo_bin_vec);
stripe_id_vec = [inference_results.stripe_id];
stripe_index = unique(stripe_id_vec);
% run some quick sanity checks on the grouping variables
mf_vec = [];
stripe_vec = [];
fluo_vec = [];
for i = 1:numel(inference_results)
    fluo_data = inference_results(i).traces;
    stripe = inference_results(i).stripe_id;
    f_bin = inference_results(i).fluo_bin;
    for j = 1:numel(fluo_data)
        mf_vec = [mf_vec nanmean(fluo_data{j})];
        stripe_vec = [stripe_vec stripe];
        fluo_vec = [fluo_vec f_bin];
    end
end

cm = jet(128);
for s = 1:numel(stripe_index)
    stripe_fig = figure;
    hold on
    for f = 1:numel(fluo_index)
        scatter(fluo_vec(fluo_vec==fluo_index(f)&stripe_vec==stripe_index(s)),mf_vec(fluo_vec==fluo_index(f)&stripe_vec==stripe_index(s)),...
            'MarkerFaceColor',cm((f-1)*10+1,:),'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
    end
    grid on
    xlabel('fluorescence bin')
    ylabel('mean trace fluorescence')
    title(['Data Integrity Check for Stripe ' num2str(stripe_index(s))])
    saveas(stripe_fig, [FigPath 'fluo_check_stripe' num2str(stripe_index(s)) '.png'])
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% iterate through results and generate summary stats %%%%%%%%%%%
hmm_results = struct;
transfer_vars = {'r','noise','soft_struct','stripe_id','fluo_bin'...
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
    r_error = r_error + 1*(min(r_sorted) > -.05*max(r_sorted));
    
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
%% generate summary table to write to csv
r_eff_vec = [hmm_results.r_eff];
kon_vec = [hmm_results.kon];
koff_vec = [hmm_results.koff];
r_mat = reshape([hmm_results.r]',numel(kon_vec),K);
ft = ~isnan(kon_vec);


inference_id_vec = 1:numel(hmm_results);
results_table = array2table([inference_id_vec' fluo_bin_vec' stripe_id_vec' kon_vec' koff_vec' ...
        r_mat r_eff_vec'], 'VariableNames',{'inf_id', 'fluo_bin', 'stripe_id',...
        'kon','koff','r1','r2','r3','r_eff'});
% save
writetable(results_table,[WritePath 'fluo_bin_results_full_final.csv'])
save([WritePath 'fluo_bin_results_full_final.mat'],'hmm_results')

%%% Generate and save list of particle IDs corresponding to each inference
particle_id_cell = cell(1,numel(inference_results));
for i = 1:numel(inference_results)
    particle_id_cell{i} = inference_results(i).particle_ids;
end
save([WritePath 'fluo_bin_particle_ids.mat'],'particle_id_cell')


% generate vector of mean fluorescent values for fluo bins
mean_fluo_vec = NaN(size(fluo_index));
for  f = 1:numel(fluo_index)
    mean_fluo_vec(f) = mean([fluo_inf_struct(tr_fluo_id==f).fluo_interp]);
end
% make condensed summary
fluo_index = unique(fluo_bin_vec);
stripe_index = unique(stripe_id_vec);
iter = 1;
for f = 1:numel(fluo_index)
    for s = 1:numel(stripe_index)
        ft = fluo_bin_vec==fluo_index(f)&stripe_id_vec==stripe_index(s)&~isnan(kon_vec);   
        ft_raw = fluo_bin_vec==fluo_index(f)&stripe_id_vec==stripe_index(s);
        % calculate statistics 
        r_med = 60*nanmean(r_eff_vec(ft));
        r_ste = 60*nanstd(r_eff_vec(ft));
        
        kon_med = 60*nanmean(kon_vec(ft));
        kon_ste = 60*nanstd(kon_vec(ft));
        
        koff_med = 60*nanmean(koff_vec(ft));
        koff_ste = 60*nanstd(koff_vec(ft));
        
        if iter == 1
            summary_mat = [fluo_index(f) mean_fluo_vec(f) stripe_index(s) sum(ft_raw) sum(ft) kon_med kon_ste koff_med koff_ste r_med r_ste];
        else
            summary_mat = vertcat(summary_mat,[fluo_index(f) mean_fluo_vec(f) stripe_index(s) sum(ft_raw) sum(ft) kon_med kon_ste koff_med koff_ste r_med r_ste]);
        end
        iter = iter + 1;
    end
end

summary_table = array2table(summary_mat(:,[1:3 6:11]), 'VariableNames',{'fluo_bin', 'mean_fluo', 'stripe_id',...
        'kon','kon_err','koff','koff_err','r','r_err'});
writetable(summary_table,[WritePath 'fluo_hmm_results_final.csv'])
%%
% plot_vars = {'kon_vec','koff_vec','cycle_vec','kon_vec./koff_vec','r_mat(:,2)','r_mat(:,3)'};
% pv_names = {'kon (s^{-1})','koff (s^{-1})','cycle time (s)', 'kon-koff ratio',...
%     'middle initiation rate (AU s^{-1})','high initiation rate (AU s^{-1})'};
% symbol_cell = {'d','o','o','o','o','o','o','o'};
% for p = 1:numel(plot_vars)
%     pv = plot_vars{p};
%     pv_num = eval(pv);
%     cmap = cm(stripe_index*15+1,:);
%     
%     fig = figure;
%     colormap(cmap);
%     hold on
%     % plot individual stripe averages first
%     lgd = {};
%     sct = [];
%     for s = stripe_index(2:end)          
%         str = ['stripe ' num2str(s)];
%         lgd = [lgd{:} {str}];
%         for f = fluo_index                    
%             sc = scatter(mean_fluo_vec(f),nanmean(pv_num(fluo_bin_vec==f&stripe_id_vec==s)),40,...
%                 s+1,'filled',symbol_cell{s+1},'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.8,'MarkerEdgeColor','black');        
%             if f == 1
%                 sct = [sct sc];
%             end
%         end
%     end
%     % now plot cross-stripe averages
%     for f = fluo_index
%         errorbar(mean_fluo_vec(f),nanmean(pv_num(fluo_bin_vec==f&stripe_id_vec==0)),...
%             nanstd(pv_num(fluo_bin_vec==f&stripe_id_vec==0)),'Color','black');
%         sc = scatter(mean_fluo_vec(f),nanmean(pv_num(fluo_bin_vec==f&stripe_id_vec==0)),...
%             40,'d','MarkerFaceColor','black','MarkerEdgeColor','black');
%         if f == 1
%             sct = [sct sc'];
%         end
%     end    
%     lgd = [lgd{:} {'all stripes'}];
%     legend(sct,lgd,'Location','best')
%     xlabel('fluorescence (AU)')
%     ylabel(pv_names{p})
%     grid on
%     title([pv_names{p} ' as a Function of Trace Fluorescence'])
% %     error('asfa')
%     % colorbar
%     saveas(fig,[FigPath 'scatter_' pv_names{p} '.png'])    
% end