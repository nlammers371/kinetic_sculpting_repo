% Script to Conduct HMM Inference on Experimental Data
close all
clear 
savio = 0; % Specify whether inference is being conducted on Savio Cluster
ap_ref_index = round((2:22)/3,1); % only used for savio inference

if savio
    addpath('/global/home/users/nlammers/repos/hmmm/src/utilities/');
    %Get environment variable from job script
    savio_groups = {str2num(getenv('SLURM_ARRAY_TASK_ID'))};    
    bin_groups = cell(1,length(savio_groups));
    for i =1:length(bin_groups)
        bin_groups{i} = ap_ref_index(savio_groups{i});
    end
else
    addpath('E:\Nick\projects\hmmm\src\utilities'); % Route to hmmm utilities folder
    bin_groups = {};
end
%-------------------------------System Vars-------------------------------%
% Core parameters
inference_times = 40*60;%(10:5:45)*60; % time at which to conduct sliding window inference
t_window = 30*60; % determines width of sliding window
K = 3; % State(s) to use for inference
w = 7; % Memory
dp_bootstrap = 1; % if 1 use bootstrap resampling at level of data points
n_bootstrap = 10; % number of bootstraps (overridden for set bootstrapping)
sample_size = 3000; % number of data points to use
min_dp_per_inf = 1000; % inference will be aborted if fewer present
project = 'eve7stripes_inf_2018_04_28';

%%%% Stable Params (these rarely change) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
Tres = 20; % Time Resolution
stop_time_inf = 60; % Specify cut-off time for inference
min_dp = 10; % min length of traces to include
clipped_ends = 1; % if one, remove final w time steps from traces
dynamic_bins = 1; % if 1, use time-resolved region classifications
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
fluo_field = 1; % specify which fluo field to (1 or 3)
% inference params
n_localEM = 25; % set num local runs
n_steps_max = 500; % set max steps per inference
eps = 1e-4; % set convergence criteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------Define Inference Variables------------------------------%

t_inf_str = num2str(round(inference_times/60));

% max num workers
if savio
    MaxWorkers = 24;
else
    MaxWorkers = 25;
end

%----------------------------Bootstrap Vars-------------------------------%
if set_bootstrap && dp_bootstrap
    error('Both set and dp bootstrapping set to 1. Choose one.')
end

if set_bootstrap
    out_string = '_set';    
elseif dp_bootstrap
    out_string = '_dp';
else
    error('inference type unspecified')
end


%----------------------------Set Write Paths------------------------------%
datapath = ['../../dat/' project '/']; %Path to raw data
% generate read and write names
dataname = ['inference_traces_' project '_dT' num2str(Tres) '.mat'];
% Load data for inference into struct named: trace_struct_final
load([datapath dataname]);

alpha = trace_struct_final(1).alpha_frac*w; % Rise Time for MS2 Loops

% Set write path (inference results are now written to external directory)
out_suffix =  ['/' project '/w' num2str(w) '_t' num2str(Tres)...
    '_alpha' num2str(round(alpha*10)) '_f' num2str(fluo_field) '_cl' num2str(clipped) ...
    '_no_ends' num2str(clipped_ends) '_tbins' num2str(dynamic_bins) ...
    '/K' num2str(K(1)) '_t_window' num2str(round(t_window/60)) '_t_inf' t_inf_str out_string '/']; 
if savio
    out_prefix = '/global/scratch/nlammers/eve7stripes_data/';
else    
    out_prefix = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
%     out_prefix = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
end
out_dir = [out_prefix out_suffix];
mkdir(out_dir);

% apply time filtering 
trace_struct_filtered = [];
for i = 1:length(trace_struct_final)
    temp_old = trace_struct_final(i);
    temp_new = struct;
    if clipped
        time = temp_old.time_interp;        
        fluo = temp_old.fluo_interp;        
    else
        time = temp_old.time_full;
        fluo = temp_old.fluo_full;                    
    end
    tf = time((time>w*Tres)&(time<stop_time_inf*60));
    if length(tf) >= min_dp
        temp_new.fluo_inf = fluo(ismember(time,tf));
        temp_new.time_inf = tf;
        temp_new.setID = temp_old.setID;
        temp_new.inference_flag = temp_old.inference_flag;
        temp_new.stripe_id_vec = temp_old.stripe_id_vec;
        temp_new.fluo_orig = temp_old.fluo;
        temp_new.time_orig = temp_old.time;
        temp_new.ParticleID = temp_old.ParticleID;
        trace_struct_filtered = [trace_struct_filtered temp_new];
    end
end
trace_struct_filtered = trace_struct_filtered([trace_struct_filtered.inference_flag]==1);  
if clipped_ends
    rm_index = 1:length(trace_struct_filtered);
    rm_indices = [];
    for i = 1:length(trace_struct_filtered)
        ft = trace_struct_filtered(i).fluo_inf;
        tt = trace_struct_filtered(i).time_inf;
        ft = ft(1:end-w-1);
        tt = tt(1:end-w-1);
        trace_struct_filtered(i).fluo_inf = ft;
        trace_struct_filtered(i).time_inf = tt;
        if length(trace_struct_filtered(i).fluo_inf) < 5
            rm_indices = [rm_indices i];
        end
    end
    trace_struct_filtered =  trace_struct_filtered(~ismember(rm_index,rm_indices));
end
set_index = [trace_struct_filtered.setID];
set_vec = unique(set_index);
if set_bootstrap % overrides n_bootstrap entry above
    n_bootstrap = length(set_vec);
end
first_time_vec = [];
last_time_vec = [];
for i = 1:length(trace_struct_filtered)
    first_time_vec = [first_time_vec min(trace_struct_filtered(i).time_inf)];
    last_time_vec = [last_time_vec max(trace_struct_filtered(i).time_inf)];
end
%% Conduct Inference
% structure array to store the analysis data
outputs = struct; % Compile outputs
local_meta = struct; % Store local_EM results
init_meta = struct; % Store initiation info

for g = 1:length(bin_groups) % loop through different AP groups
    bin_list = bin_groups{g}; % get groups for present iteration                        
    for t = 1:length(inference_times)
        t_inf = inference_times(t);
        t_start = t_inf - t_window/2;
        t_stop = t_inf + t_window/2;
        for b = 1:n_bootstrap
            iter_start = now;
            s = (g-1)*n_bootstrap + b;
            local_struct = struct;
            init_struct = struct;
            output = struct;
            % if using set bootstrapping, select set to hold out
            if set_bootstrap
                boot_set = set_vec(b);
            end
            % Use current time as unique inference identifier 
            inference_id = num2str(round(10e5*now));
            if set_bootstrap
                inference_id = [inference_id '_set' num2str(boot_set)];
            end
            % Generate filenames            
            fName_sub = ['eveSet_w' num2str(w) '_K' num2str(K) ...
                '_bin' num2str(round(10*bin_list(1))) '_' num2str(round(10*bin_list(end))) ...
                '_time' num2str(round(t_inf/60)) '_t' inference_id];                
            out_file = [out_dir '/' fName_sub];
            
            % Extract fluo_data
            temp_filter = (last_time_vec-5*60) > t_start ...
                    & (first_time_vec+5*60) <= t_stop;
            stripe_id_vec = NaN(1,length(temp_filter));
            for m = 1:length(trace_struct_filtered)
                if temp_filter(m)==1
                    tt = trace_struct_filtered(m).time_orig;
                    ft = trace_struct_filtered(m).fluo_orig;
                    tt = tt(~isnan(ft));
                    time_filter = tt>=t_start & tt < t_stop;
                    stripe_id_trace = trace_struct_filtered(m).stripe_id_vec(time_filter);
                    stripe_id_trace(isnan(stripe_id_trace)) = Inf;
                    stripe_id_vec(m) = round(mode(trace_struct_filtered(m).stripe_id_vec(time_filter)),1);
                end
            end                
            trace_filter = ismember(stripe_id_vec,bin_list);
            if set_bootstrap
                trace_ind = find(([trace_struct_filtered.setID]~=boot_set)&...
                    trace_filter);
            else
                trace_ind = find(trace_filter);
            end
            % generate structure containing only elligible trace
            % fragments
            inference_set = [];            
            for m = 1:length(trace_ind)
                temp = trace_struct_filtered(trace_ind(m));
                tt = temp.time_inf;
                ft = temp.fluo_inf;
                temp.time_inf = tt(tt>=t_start & tt < t_stop);
                temp.fluo_inf = ft(tt>=t_start & tt < t_stop);
                if sum(temp.fluo_inf>0) > 1 % remove pure zero fragments
                    inference_set = [inference_set temp];
                end                
            end            
            if isempty(inference_set)
                skip_flag = 1;
            else                
                set_size = length([inference_set.time_inf]);                                
                if set_size < min_dp_per_inf
                    skip_flag = 1;
                else
                    skip_flag = 0;
                end
            end
            if skip_flag
                warning('Too few data points. Skipping')                
            else 
                sample_index = 1:length(inference_set);
                particle_id_list = [];
                particle_time_cell = {};
                if dp_bootstrap                        
                    ndp = 0;    
                    sample_ids = [];                    
                    %Reset bootstrap size to be on order of set size for small bins
                    if set_size < sample_size
                        sample_size = ceil(set_size/1000)*1000;
                    end
                    while ndp < sample_size
                        tr_id = randsample(sample_index,1);
                        sample_ids = [sample_ids tr_id];
                        particle_id_list = [particle_id_list inference_set(tr_id).ParticleID];
                        particle_time_cell = [particle_time_cell{:} {inference_set(tr_id).time_inf}];
                        ndp = ndp + length(inference_set(tr_id).time_inf);
                    end
                    fluo_data = cell([length(sample_ids), 1]);                
                    sample_particles = [inference_set(sample_ids).ParticleID];
                    for tr = 1:length(sample_ids)
                        fluo_data{tr} = inference_set(sample_ids(tr)).fluo_inf;                    
                    end                                
                else % Take all relevant traces if not bootstrapping
                    error('non-bootstrap option is deprecated')
                end                 
                % Random initialization of model parameters
                param_init = initialize_random (K, w, fluo_data);
                % Approximate inference assuming iid data for param initialization                
                local_iid_out = local_em_iid_reduced_memory_truncated (fluo_data, param_init.v, ...
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
                output.total_time = 100000*(now - iter_start);            
                % Save inference ID variables
                output.stripe_id = bin_list;
                output.boot_set = NaN;
                if set_bootstrap
                    output.boot_set = boot_set;
                end             
                output.t_window = t_window;
                output.t_inf = t_inf;
                output.fluo_type = fluo_field;
                output.dp_bootstrap_flag = dp_bootstrap;
                output.set_bootstrap_flag = set_bootstrap;
                output.iter_id = b;
                output.start_time_inf = 0;
                output.stop_time_inf = stop_time_inf;                            
                output.clipped = clipped;
                if dp_bootstrap || set_bootstrap
%                         output.traces = sample_ids;    
                    output.particle_ids = particle_id_list;
                    output.particle_times = particle_time_cell;
                    output.traces = fluo_data;
                    output.N = ndp;
                end
                % Other inference characteristics            
                output.w = w;
                output.alpha = alpha;
                output.deltaT = Tres;  
                output.dynamic_bins = 1;
            end
            output.skip_flag = skip_flag;
            output.inference_traces = fluo_data;
            save([out_file '.mat'], 'output');           
        end  
    end
end    
