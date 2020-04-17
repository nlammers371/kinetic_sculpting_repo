% Script to Conduct HMM Inference on Experimental Data
close all
clear 
savio = 1; % Specify whether inference is being conducted on Savio Cluster
stripe_ref = 0:7; % only used for savio inference

if savio
    addpath('/global/home/users/nlammers/repos/hmmm/src/utilities/');
    %Get environment variable from job script
    savio_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};    
    stripe_id_index = NaN(1,length(savio_groups));
    for i =1:length(stripe_id_index)
        stripe_id_index(i) = stripe_ref(savio_groups{i});
    end
else
    addpath('E:\Nick\projects\hmmm\src\utilities'); % Route to hmmm utilities folder
    stripe_id_index = stripe_ref;
end

%-------------------------------System Vars-------------------------------%
% Core parameters
K = 3; % State(s) to use for inference
w = 7; % Memory
dp_bootstrap = 1; % if 1 use bootstrap resampling at level of data points
n_bootstrap = 10; % number of bootstraps (overridden for set bootstrapping)
sample_size = 5000; % number of data points to use
min_dp_per_inf = 1000; % inference will be aborted if fewer present
project = 'revision_fluo_bins';
ReadPath = '../../dat/revisions/';
use_inf_tape = 1;
% check to see if inference tape structure exists
if isfile(['../../out/revisions/' project '/inference_tape.mat']) && use_inf_tape
    load(['../../out/revisions/' project '/inference_tape.mat']);
    n_boot_vec = flipud(inference_tape.n_boot_vec);
    stripe_id_vec = flipud(inference_tape.stripe_id_vec);
    fluo_bin_vec = flipud(inference_tape.fluo_bin_vec);
else
    f_bins = 1:10;
    fluo_bin_vec = repelem(f_bins,numel(stripe_id_index));
    stripe_id_vec = repelem(stripe_id_index,numel(f_bins));
    n_boot_vec = repelem(n_bootstrap,numel(stripe_id_vec));
end
    
%%%% Stable Params (these rarely change) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
Tres = 20; % Time Resolution
min_dp = 10; % min length of traces to include

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
load([ReadPath 'fluo_inf_struct.mat']);
alpha = fluo_inf_struct(1).alpha_frac*w; % Rise Time for MS2 Loops

% Set write path (inference results are now written to external directory)
out_suffix =  ['/' project '/w' num2str(w) '_K' num2str(K) '_' project '/']; 
if savio
    out_prefix = '/global/scratch/nlammers/eve7stripes_data/';
else    
    out_prefix = '../../out/revisions/';
%     out_prefix = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
end
out_dir = [out_prefix out_suffix];
mkdir(out_dir);


%% Conduct Inference
% iterate through inference groups
for f = 1:numel(fluo_bin_vec)         
    stripe_bin = stripe_id_vec(f); % get groups for present iteration   
    n_bootstrap = n_boot_vec(f);
    for b = 1:n_bootstrap                
        local_struct = struct;
        init_struct = struct;
        output = struct;        
        % Use current time as unique inference identifier 
        inference_id = num2str(round(10e5*now));        
        % Generate filenames            
        fName_sub = ['eveSet_w' num2str(w) '_K' num2str(K) ...
            '_stripe' num2str(stripe_bin) '_fbin' num2str(fluo_bin_vec(f)) '_t' inference_id];                
        out_file = [out_dir '/' fName_sub];
        % Extract fluo_data        
        stripe_ids = round([fluo_inf_struct.stripe_id]);
        fluo_bins = [fluo_inf_struct.FluoBin];
        if stripe_bin ~= 0
            trace_filter = ismember(stripe_ids,stripe_bin)&fluo_bins==fluo_bin_vec(f);
        else
            trace_filter = fluo_bins==fluo_bin_vec(f);
        end
        trace_ind = find(trace_filter);        
        % generate structure containing only elligible trace        
        inference_set = fluo_inf_struct(trace_ind);                                    
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
                sample_particles = [inference_set(sample_ids).ParticleID];
                for tr = 1:length(sample_ids)
                    fluo_data{tr} = inference_set(sample_ids(tr)).fluo_interp;                    
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
            output.stripe_id = stripe_bin;                                
            output.fluo_bin = fluo_bin_vec(f);                                                
            output.dp_bootstrap_flag = dp_bootstrap;               
            output.iter_id = b;
            output.traces = fluo_data;
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
%             output.inference_traces = fluo_data;
        save([out_file '.mat'], 'output');           
    end      
end    
