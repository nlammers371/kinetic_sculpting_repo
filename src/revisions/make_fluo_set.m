% script to generate fluorescence-dependent binning for inference
clear 
close all

OutPath = '../../dat/revisions/';
% generate read and write names
DataName = 'inference_traces_w_stripe_ids.mat';
% Load data for inference into struct named: trace_struct_final
load([OutPath DataName]);

% bin traces by mean fluorescence
n_bins = 10;
min_obs = 10;
mf_vec = NaN(1,numel(trace_struct_final));
% generate indexing vectors for grouping
pt_vec = [trace_struct_final.ParticleID];
stripe_vec = NaN(size(pt_vec));
for i = 1:numel(trace_struct_final)
    f_vec = trace_struct_final(i).fluo_interp;
    s_id_vec = trace_struct_final(i).mike_stripe_id;
    if numel(f_vec) >= min_obs && max(isnan(s_id_vec))==0
        mf_vec(i) = mean(f_vec);
        stripe_vec(i) = round(mean(s_id_vec));
    end
end
% calculate activity percentiles and group traces accordingly
prctile_vec = linspace(100/n_bins, 100,n_bins);
f_id_vec_full = NaN(size(stripe_vec));
fluo_pct_vec =  NaN(size(prctile_vec));
for i = 1:numel(prctile_vec)
    pct = prctile(mf_vec,prctile_vec(i));
    f_id_vec_full(mf_vec<=pct&isnan(f_id_vec_full)) = i;
    fluo_pct_vec(i) = pct;
end
% generate trace structure
i_pass = 1;
fluo_inf_struct = struct;
for i = 1:numel(pt_vec)    
    if ~isnan(mf_vec(i))
        fluo_inf_struct(i_pass).ParticleID = pt_vec(i);
        fluo_inf_struct(i_pass).stripe_id = stripe_vec(i);
        fluo_inf_struct(i_pass).FluoBin = f_id_vec_full(i);
        fluo_inf_struct(i_pass).fluo_interp = trace_struct_final(i).fluo_interp;
        fluo_inf_struct(i_pass).alpha_frac = trace_struct_final(i).alpha_frac;    
        i_pass = i_pass + 1;
    end
end
save([OutPath 'fluo_inf_struct.mat'] , 'fluo_inf_struct')