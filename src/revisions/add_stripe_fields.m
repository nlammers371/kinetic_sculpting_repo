% script to incorporate Mike's strip classification scheme 
clear 
close all
% basic ID params
K = 3; % State(s) to use for inference
w = 7; % Memory
project = 'eve7stripes_inf_2018_04_28'; % project identifier
Tres = 20; % Time Resolution
% data paths
datapath = ['../../dat/' project '/']; %Path to raw data
OutPath = '../../dat/revisions/';
mkdir(OutPath)
% generate read and write names
dataname = ['inference_traces_' project '_dT' num2str(Tres) '.mat'];
% load data for inference into struct named: trace_struct_final
load([datapath dataname]);
% load Mike's stripe assigment key
stripe_assignments = readtable('E:\Nick\projects\all_the_stripes_cd\dat\revisions\stripe_assignments_eve_9jan19.csv');

add_particle_id_vec = stripe_assignments.particle_id;
% add stripe info to data structure
tr_particle_id_vec = [trace_struct_final.ParticleID];
for i = 1:numel(tr_particle_id_vec)
    ParticleID = tr_particle_id_vec(i);
    old_stripe = trace_struct_final(tr_particle_id_vec==ParticleID).stripe_id_vec_interp;
    ft = add_particle_id_vec==ParticleID;
    if sum(ft) ~= 0
        new_stripe = stripe_assignments(ft,:).stripe_id;
        stripe_distance = stripe_assignments(add_particle_id_vec==ParticleID,:).burst_stripe_ap_d;
        if numel(old_stripe)~=numel(new_stripe)        
            error('adsa')
        end
    else
        new_stripe = NaN(size(old_stripe));
        stripe_distance = NaN(size(old_stripe));
    end
    trace_struct_final(tr_particle_id_vec==ParticleID).burst_stripe_ap_d = stripe_distance';
    trace_struct_final(tr_particle_id_vec==ParticleID).mike_stripe_id = new_stripe';
end

save([OutPath 'inference_traces_w_stripe_ids.mat'],'trace_struct_final')