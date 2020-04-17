% Script to Conduct HMM Inference on Experimental Data
close all
clear 

%-------------------------------System Vars-------------------------------%
% Core parameters
K = 3; % State(s) to use for inference
w = 7; % Memory
dp_bootstrap = 1; % if 1 use bootstrap resampling at level of data points
n_bootstrap = 5; % number of bootstraps (overridden for set bootstrapping)
sample_size = 5000; % number of data points to use
min_dp_per_inf = 1000; % inference will be aborted if fewer present
project = 'revision_fluo_bins_v2';
ReadPath = '../../dat/revisions/';
savio=1; % Specify whether inference is being conducted on Savio Cluster
stripe_id_index = 0:7; % only used for savio inference
t_start = 0*60; % minimum time for inclusion in inference

% add path to utilities folder
if savio
    addpath('../utilities/');
    %Get environment variable from job script
else
    addpath('../utilities'); % Route to hmmm utilities folder
end
    
%%%% Stable Params (these rarely change) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
Tres = 20; % Time Resolution
minDP = 10; % min length of traces to include

% inference params
n_localEM = 25; % set num local runs
n_steps_max = 500; % set max steps per inference
eps = 1e-4; % set convergence criteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------Define Inference Variables------------------------------
% max num workers
if savio
    MaxWorkers = 24;
else
    MaxWorkers = 25;
end

%----------------------------Set Write Paths------------------------------%
load([ReadPath 'inference_traces_w_stripe_ids.mat']);
alpha = trace_struct_final(1).alpha_frac*w; % Rise Time for MS2 Loops

% Set write path (inference results are now written to external directory)
out_suffix =  ['/' project '/w' num2str(w) '_K' num2str(K) '_t' num2str(t_start/60) '_' project '/']; 
if savio
    out_prefix = '/global/scratch/nlammers/eve7stripes_data_v2/';
else    
    out_prefix = '../../out/revisions/';
%     out_prefix = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
end
out_dir = [out_prefix out_suffix];
mkdir(out_dir);

% filter traces and make fluorescence bins
trace_struct_filtered = [];
for i = 1:length(trace_struct_final)
    temp = struct;
    time = trace_struct_final(i).time_interp(w+1:end);    
    fluo = trace_struct_final(i).fluo_interp(w+1:end);  
    time_ft = time  >= t_start;    
    if sum(time_ft) >= minDP
        temp.fluo_interp = fluo(time_ft);
        temp.time_interp = time(time_ft);             
        temp.Stripe = mean(round(trace_struct_final(i).mike_stripe_id(time_ft)));           
        temp.MeanFluo = nanmean(fluo(time_ft));
        temp.ParticleID = trace_struct_final(i).ParticleID;    
        temp.N = sum(time_ft);
        trace_struct_filtered = [trace_struct_filtered temp];    
    end
end

% generate fluorescence bins that are consistent across all stripes
f_bins = 1:9;
stripe_id_vec = [trace_struct_filtered.Stripe];
mean_fluo_vec = [trace_struct_filtered.MeanFluo];

nTotal = sum([trace_struct_filtered.N]);
prctile_vec = linspace(0,1,numel(f_bins)+1);
% get quantile bins
fluo_quantiles = quantile(mean_fluo_vec,prctile_vec);
fluo_id_vec = discretize(mean_fluo_vec,fluo_quantiles);
for i = 1:numel(trace_struct_filtered)
    trace_struct_filtered(i).FluoBin = fluo_id_vec(i);
end
stripe_vec_inf = repelem(stripe_id_index,numel(f_bins));
fluo_vec_inf = repelem(f_bins,numel(stripe_id_index));


% Conduct Inference
% iterate through inference groups
rng('shuffle')
inference_list = randsample(1:numel(fluo_vec_inf),numel(fluo_vec_inf),false);
for f = inference_list         
    stripe_bin = stripe_vec_inf(f); % get groups for present iteration      
    for b = 1:n_bootstrap                
        local_struct = struct;
        init_struct = struct;
        output = struct;        
        
        % Use current time as unique inference identifier 
        inference_id = num2str(round(10e5*now));        
        
        % Generate filenames            
        fName_sub = ['eveSet_w' num2str(w) '_K' num2str(K) ...
            '_stripe' num2str(stripe_bin) '_fbin' num2str(fluo_vec_inf(f)) '_t' inference_id];                
        out_file = [out_dir '/' fName_sub];
        
        % Extract fluo_data        
        stripe_ids = round([trace_struct_filtered.Stripe]);
        fluo_bins = [trace_struct_filtered.FluoBin];
        if stripe_bin ~= 0
            trace_filter = ismember(stripe_ids,stripe_bin)&fluo_bins==fluo_vec_inf(f);
        else
            trace_filter = fluo_bins==fluo_vec_inf(f);
        end
        trace_ind = find(trace_filter);        
        
        % generate structure containing only elligible trace        
        inference_set = trace_struct_filtered(trace_ind);                                    
        set_size = length([inference_set.fluo_interp]);   
        
        % check that we have enough data points
        skip_flag = set_size < min_dp_per_inf;

        if skip_flag
            warning('Too few data points. Skipping')                
        else 
            sample_index = 1:length(inference_set);
            particle_id_list = [];
            if dp_bootstrap                        
                ndp = 0;    
                sample_ids = [];                    
                %Reset bootstrap size to be on order of set size for small bins
                samp_size = sample_size;
                if set_size < samp_size
                    samp_size = ceil(set_size/1000)*1000;
                end
                while ndp < samp_size
                    tr_id = randsample(sample_index,1);
                    sample_ids = [sample_ids tr_id];
                    particle_id_list = [particle_id_list inference_set(tr_id).ParticleID];                        
                    ndp = ndp + length(inference_set(tr_id).fluo_interp);
                end
                fluo_data = cell([length(sample_ids), 1]);  
                time_data = cell([length(sample_ids), 1]);  
                sample_particles = [inference_set(sample_ids).ParticleID];
                for tr = 1:length(sample_ids)
                    fluo_data{tr} = inference_set(sample_ids(tr)).fluo_interp;   
                    time_data{tr} = inference_set(sample_ids(tr)).time_interp;   
                end                                
            else % Take all relevant traces if not bootstrapping
                error('non-bootstrap option is deprecated')
            end                 
            % Random initialization of model parameters
            param_init = initialize_random (K, w, fluo_data);
            % Approximate inference assuming iid data for param initialization                
            local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                param_init.noise, K, w, alpha, n_steps_max, eps);
            noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
            v_iid = exp(local_iid_out.v_logs);            
            p = gcp('nocreate');
            if isempty(p)
                parpool(MaxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
            elseif p.NumWorkers > MaxWorkers
                delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                parpool(MaxWorkers);
            end
            parfor i_local = 1:n_localEM % Parallel Local EM                
                % Random initialization of model parameters
                param_init = initialize_random_with_priors(K, noise_iid, v_iid);
                % Get Intial Values
                pi0_log_init = log(param_init.pi0);
                A_log_init = log(param_init.A);
                v_init = param_init.v;
                noise_init = param_init.noise;
                % Record
                init_struct(i_local).A_init = exp(A_log_init);                
                init_struct(i_local).v_init = v_init;
                init_struct(i_local).noise_init = noise_init;                
                init_struct(i_local).subset_id = i_local;
                %--------------------LocalEM Call-------------------------%
                local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                    v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                    alpha, n_steps_max, eps);                    
                %---------------------------------------------------------%                
                % Save Results 
                local_struct(i_local).inference_id = inference_id;
                local_struct(i_local).subset_id = i_local;
                local_struct(i_local).logL = local_out.logL;                
                local_struct(i_local).A = exp(local_out.A_log);
                local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / Tres;                                
                local_struct(i_local).noise = 1/exp(local_out.lambda_log);
                local_struct(i_local).pi0 = exp(local_out.pi0_log);
                local_struct(i_local).total_steps = local_out.n_iter;               
                local_struct(i_local).soft_struct = local_out.soft_struct;               
            end
            [logL, max_index] = max([local_struct.logL]); % Get index of best result           
            % Save parameters from most likely local run
            output.pi0 =local_struct(max_index).pi0;                        
            output.r = local_struct(max_index).r(:);
            output.noise = local_struct(max_index).noise;
            output.A = local_struct(max_index).A(:);
            output.A_mat = local_struct(max_index).A;            
            % get soft-decoded structure
            output.soft_struct = local_struct(max_index).soft_struct;
            % Info about run time
            output.total_steps = local_struct(max_index).total_steps;                                                  
            % Save inference ID variables
            output.Stripe = stripe_bin;                                
            output.fluo_bin = fluo_vec_inf(f);                                                
            output.dp_bootstrap_flag = dp_bootstrap;               
            output.iter_id = b;
            output.fluo_traces = fluo_data;
            output.time_traces = time_data;
            output.particle_ids = particle_id_list;
            if dp_bootstrap 
                output.N = ndp;
            end
            % Other inference characteristics            
            output.w = w;
            output.alpha = alpha;
            output.deltaT = Tres;                  
        end
        output.skip_flag = skip_flag;
        save([out_file '.mat'], 'output');           
    end      
end    
