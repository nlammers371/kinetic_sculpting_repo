% script to assess misssing pieces of fluo bin inference
clear
close all
% id variables 
project = 'revision_fluo_bins';
readPath = ['../../out/revisions/' project '/'];
load('../../dat/revisions/fluo_inf_struct.mat');
summary_table = readtable([readPath 'fluo_bin_summary.csv']);
% analysis params
n_boots_min = 10;
n_boot_vec = summary_table.n_dup_final;
fluo_bin_vec = summary_table.fluo_bin;
stripe_id_vec = summary_table.stripe_id;
kon_vec = summary_table.kon_med;
% generate list of stripe-fluo combinations that need additional data
% points
deficient_ids = find(n_boot_vec < n_boots_min & ~isnan(kon_vec));
inference_tape = struct;
inference_tape.n_boot_vec = n_boots_min - n_boot_vec(deficient_ids);
inference_tape.stripe_id_vec = stripe_id_vec(deficient_ids);
inference_tape.fluo_bin_vec = fluo_bin_vec(deficient_ids);
% save
save([readPath 'inference_tape.mat'],'inference_tape')