% Script to call inference compilation function
%%% Define necessary variables
% Script to Conduct HMM Inference on Experimental Data
close all
clear 
addpath('../utilities');
%-------------------------------System Vars-------------------------------%
w = 7; % Memory
K = 3; % State(s) to use for inference
clipped_ends = 1; % if one, remove final w time steps from traces
dynamic_bins = 1; % if 1, use time-resolved region classifications
clipped = 1; % if 0 use "full" trace with leading and trailing 0's
fluo_type = 1; % specify which fluo field to (1 or 3)
t_window = 30; % determines width of sliding window
t_inf = 40;
inference_type = 'dp';
project = 'eve7stripes_inf_2018_04_28';
Tres = 20;
alpha = 1.4;
% DPFolder = 'E:/Nick/Dropbox (Garcia Lab)/eve7stripes_data/inference_out/';
DPFolder = 'S:\Nick\GarciaLab\Dropbox (Garcia Lab)\eve7stripes_data\inference_out\';
%-----------------------------ID Variables--------------------------------%
stripe_range = 1:7;
bin_range_vec = [];
for i = 1:length(stripe_range)
    for j = 1:3
        bin_range_vec = stripe_range(i) - j/3 + 2/3;
    end
end
%Generate filenames and writepath
id_var = [ '/w' num2str(w) '_t' num2str(Tres) '_alpha' num2str(round(alpha*10)) ...
    '_f' num2str(fluo_type) '_cl' num2str(clipped) '_no_ends' num2str(clipped_ends) ...
    '_tbins' num2str(dynamic_bins) '/K' num2str(K) '_t_window' num2str(t_window) ...
     '_t_inf' num2str(t_inf) '_' inference_type '/']; 

folder_path =  [DPFolder  project  id_var ];
OutPath = ['../../dat/' project  id_var];
FigPath = ['../../fig/experimental_system/' project '/' id_var];
mkdir(OutPath)
mkdir(FigPath)

%---------------------------------Read in Files---------------------------%


file_list = dir([folder_path '\' '*.mat']);
inference_results = struct;
f_pass = 1;
for f = 1:numel(file_list)
    load([file_list(f).folder '\' file_list(f).name]);
    fn = fieldnames(output);
    if numel(fn)>3
        for i = 1:numel(fn)
            inference_results(f_pass).(fn{i}) = output.(fn{i});
        end
        f_pass = f_pass + 1
    end
end
%%
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

%% Iterate through result sets and concatenate into 1 combined struct
glb_all = struct;
f_pass = 1;
% p_mat = NaN(length(filenames),20);
% dup_list = [];
for f = 1:length(filenames)
    % load the eve validation results into a structure array 'output'    
    load([folder_path filenames{f}]);
%     if output.skip_flag == 1
    if length(fieldnames(output)) < 2 
        continue
    end    
    for fn = fieldnames(output)'
        glb_all(f_pass).(fn{1}) = output.(fn{1});
    end
    glb_all(f_pass).source = filenames{f};            
    f_pass = f_pass + 1;    
end
% dup_list = unique(dup_list);
% %%
%%% ------------------------------Fig Calculations-------------------------%
%Adjust rates as needed (address negative off-diagonal values)
%Define convenience Arrays and vectors
alpha = glb_all(1).alpha;
bin_vec = [];
bin_cell = {};
for i = 1:length(glb_all)
    if min(glb_all(i).stripe_id) == .7 && max(glb_all(i).stripe_id) == 7.3
        bin_vec = [bin_vec 0];
    else
        bin_vec = [bin_vec mean(glb_all(i).stripe_id)];
    end
    bin_cell = [bin_cell{:} {glb_all(i).stripe_id}];
end
% error('afs')
% bin_vec = [glb_all.APbin];
bin_range = unique(bin_vec);
time_vec = [glb_all.t_inf];
time_index = unique(time_vec);
initiation_rates = zeros(K,length(glb_all)); % r
pi0_all = zeros(K,length(glb_all));    
noise_all = zeros(1,length(glb_all)); % sigma  
n_dp_all = zeros(1,length(glb_all));
n_traces_all = zeros(1,length(glb_all));
% effective_times_all = zeros(1,length(glb_all));
dwell_all = NaN(K,length(glb_all));
quality_flags = ones(1,length(glb_all)); % if 1 will be used for summary stats
% Extract variables and perform rate fitting as needed
for i = 1:length(glb_all)    
    [r_inf, ranked_r] = sort(60*[glb_all(i).r]); 
    if min(r_inf) < -500 % address structural issues with inference
        quality_flags(i) = 0;       
    end
    initiation_rates(:,i) = r_inf;
    pi0 = glb_all(i).pi0(ranked_r);    
    noise_all(i) = sqrt(glb_all(i).noise); 
    n_dp_all(i) = glb_all(i).N;
    n_traces_all(i) = length(glb_all(i).particle_ids);
%     effective_times_all(i) = nanmean([glb_all(i).time_data{:}]);
    A = reshape(glb_all(i).A,K,K);
    A = A(ranked_r, ranked_r);
    %Obtain raw R matrix
    R = prob_to_rate(A,Tres);
    glb_all(i).R_mat = R;
    Rcol = reshape(R,1,[]);
    glb_all(i).R = Rcol;
    %Check for imaginary and negative elements. If present, perform rate
    %fit   
    glb_all(i).r_fit_flag = 0;    
    if ~isreal(Rcol) || sum(Rcol<0)>K
        glb_all(i).r_fit_flag = 1;
        out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
        glb_all(i).R_fit = 60*out.R_out;
        r_diag = 60*diag(out.R_out);
    elseif false %sum(Rcol<0)>K
        glb_all(i).r_fit_flag = 2;
        R_conv = 60*(R - eye(K).*R);
        R_conv(R_conv<0) = 0;
        R_conv(eye(K)==1) = -sum(R_conv);
        glb_all(i).R_fit = R_conv;
        r_diag = diag(R_conv);
    else
        glb_all(i).R_fit = 60*glb_all(i).R_mat;
        r_diag = 60*diag(glb_all(i).R_mat); 
    end    
    dwell_all(:,i) = -1./r_diag;
end

% Save rate info
%Make 3D Arrays to stack matrices
R_orig_array = NaN(length(glb_all),K^2);
R_fit_array = NaN(length(glb_all),K^2);
A_array = NaN(length(glb_all),K^2);
for i = 1:length(glb_all)
    R_fit_array(i,:) = reshape(glb_all(i).R_fit,1,[]);
    R_orig_array(i,:) = reshape(real(glb_all(i).R_mat),1,[]);
    A_array(i,:) = reshape(real(glb_all(i).A_mat),1,[]);
end

% 2D Arrays to store moments (A's are calculated to have for output structure)
med_R_orig = zeros(length(time_index),K^2,length(bin_range));
std_R_orig = zeros(length(time_index),K^2,length(bin_range));
avg_A = zeros(length(time_index),K^2,length(bin_range));
% std_A = zeros(length(time_index),K^2,length(bin_range));
med_R_fit = zeros(length(time_index),K^2,length(bin_range));
std_R_fit = zeros(length(time_index),K^2,length(bin_range));
med_on_ratio = NaN(length(time_index),length(bin_range));
std_on_ratio = NaN(length(time_index),length(bin_range));
med_off_ratio = NaN(length(time_index),length(bin_range));
std_off_ratio = NaN(length(time_index),length(bin_range));
for b = 1:length(bin_range)
    bin = bin_range(b);
    for t = 1:length(time_index)
        time = time_index(t);
        q_filter = quality_flags&bin_vec==bin&time_vec==time;
        A_mean = mean(A_array(q_filter,:),1);
        A_reshape = reshape(A_mean,K,K);
        %Normalize A
        A_mean = reshape(A_reshape ./ repmat(sum(A_reshape),K,1),1,[]);
        avg_A(t,:,b) = A_mean;    
        %Calculate and Store R moments        
        
        med_R_orig(t,:,b) = median(R_orig_array(q_filter,:),1);        
        std_R_orig(t,:,b) = .5*(quantile(R_orig_array(q_filter,:),.75,1)-...
                                quantile(R_orig_array(q_filter,:),.25,1));            
        med_R_fit(t,:,b) = median(R_fit_array(q_filter,:),1);
        std_R_fit(t,:,b) = .5*(quantile(R_fit_array(q_filter,:),.75,1)-...
                               quantile(R_fit_array(q_filter,:),.25,1));
        if K == 3
            k_on_rat = R_fit_array(q_filter,2)./...
                            R_fit_array(q_filter,6);
            med_on_ratio(t,b) = median(k_on_rat);
            std_on_ratio(t,b) = .5*(quantile(k_on_rat,.75,1)-...
                                   quantile(k_on_rat,.25,1));
            k_off_rat = R_fit_array(q_filter,8)./...
                            R_fit_array(q_filter,4);
            med_off_ratio(t,b) = median(k_off_rat);
            std_off_ratio(t,b) = .5*(quantile(k_off_rat,.75,1)-...
                                   quantile(k_off_rat,.25,1));
        end
    end
end

%Occupancy
occupancy = zeros(K,length(glb_all));
for i = 1:length(glb_all)
    [~, ranked_r] = sort([glb_all(i).r]);
    A = glb_all(i).A_mat;
    A = A(ranked_r,ranked_r);
    [V,D] = eig(A);
    ind = diag(D)==max(diag(D));
    steady = V(:,ind)./sum(V(:,ind));
    occupancy(:,i) = steady;
    glb_all(i).occupancy = steady;
end
med_dwell = NaN(K,length(bin_range),length(time_index));
std_dwell = NaN(K,length(bin_range),length(time_index));
med_initiation = NaN(K,length(bin_range),length(time_index));
std_initiation = NaN(K,length(bin_range),length(time_index));
med_init_ratio = NaN(length(time_index),length(bin_range));
std_init_ratio = NaN(length(time_index),length(bin_range));
med_occupancy = NaN(K,length(bin_range),length(time_index));
std_occupancy = NaN(K,length(bin_range),length(time_index));
avg_pi0 = NaN(K,length(bin_range),length(time_index));
std_pi0 = NaN(K,length(bin_range),length(time_index));
med_noise = NaN(1,length(bin_range),length(time_index));
std_noise = NaN(1,length(bin_range),length(time_index));
med_n_dp = NaN(1,length(bin_range),length(time_index));
med_n_tr = NaN(1,length(bin_range),length(time_index));
n_boots_used = NaN(length(time_index),length(bin_range));
n_boots_total = NaN(length(time_index),length(bin_range));
% med_eff_time = NaN(1,length(bin_range),length(time_index));
for i = 1:length(bin_range)
    bin = bin_range(i);
    for t = 1:length(time_index)
        time = time_index(t);
        q_filter = quality_flags&bin_vec==bin&time_vec==time;
        n_boots_used(t,i) = sum(q_filter);
        n_boots_total(t,i) = sum(bin_vec==bin&time_vec==time);
        for k = 1:K            
            if ~isempty(initiation_rates(k,q_filter))
                med_initiation(k,i,t) =  median(initiation_rates(k,q_filter));            
                std_initiation(k,i,t) = .5*(quantile(initiation_rates(k,q_filter),.75)-...
                        quantile(initiation_rates(k,q_filter),.25));                  
                med_occupancy(k,i,t) = median(occupancy(k,q_filter));
                std_occupancy(k,i,t) = .5*(quantile(occupancy(k,q_filter),.75)-...
                        quantile(occupancy(k,q_filter),.25)); 
                med_dwell(k,i,t) = median(dwell_all(k,q_filter));
                std_dwell(k,i,t) = .5*(quantile(dwell_all(k,q_filter),.75)-...
                        quantile(dwell_all(k,q_filter),.25));        
    %             avg_pi0(k,i,t) = mean(pi0_all(k,q_filter));
    %             std_pi0(k,i,t) = std(pi0_all(k,q_filter));
            end
        end     
        if K == 3
            init_rat = initiation_rates(3,q_filter)./...
                            initiation_rates(2,q_filter);
            med_init_ratio(t,i) = median(init_rat);
            std_init_ratio(t,i) = .5*(quantile(init_rat,.75)-...
                        quantile(init_rat,.25));
        end
        med_noise(i,t) = median(noise_all(q_filter));
        std_noise(i,t) = .5*(quantile(noise_all(q_filter),.75)-...
                    quantile(noise_all(q_filter),.25));
        med_n_dp(i,t) = median(n_dp_all(q_filter));
        med_n_tr(i,t) = median(n_traces_all(q_filter));
%         med_eff_time(i,t) = median(effective_times_all(q_filter));
    end
end
% off rate, on rate and emission rates 
if K == 3
    on_rate_eff_med = NaN(length(time_index),length(bin_range));
    on_rate_eff_ste = NaN(length(time_index),length(bin_range));
    off_rate_eff_med = NaN(length(time_index),length(bin_range));
    off_rate_eff_ste = NaN(length(time_index),length(bin_range));
    init_rate_eff_med = NaN(length(time_index),length(bin_range));
    init_rate_eff_ste = NaN(length(time_index),length(bin_range));

    for b = 1:length(bin_range)
        bin = bin_range(b);
        for t = 1:length(time_index)
            time = time_index(t);     
            q_filter = quality_flags&bin_vec==bin&time_vec==time;
            %Calculate and Store R moments  
            r21_vec = R_fit_array(q_filter,2)/2;
            r12_vec = R_fit_array(q_filter,4);
            r32_vec = R_fit_array(q_filter,6);
            r23_vec = R_fit_array(q_filter,8)/2;
            pi1_vec = occupancy(1,q_filter)';
            pi2_vec = occupancy(2,q_filter)';
            pi3_vec = occupancy(3,q_filter)';
            on_rate_all = (r21_vec.*pi1_vec+r32_vec.*(r32_vec./(r32_vec+r12_vec)).*pi2_vec)./(pi1_vec+(r32_vec./(r32_vec+r12_vec)).*pi2_vec);
            on_rate_eff_med(t,b) = median(on_rate_all);
            on_rate_eff_ste(t,b) = .5*(quantile(on_rate_all,.75)-quantile(on_rate_all,.25));
            off_rate_all = (r23_vec.*pi3_vec+r12_vec.*(r12_vec./(r32_vec+r12_vec)).*pi2_vec)./(pi3_vec+(r12_vec./(r32_vec+r12_vec)).*pi2_vec);
            off_rate_eff_med(t,b) = median(off_rate_all);
            off_rate_eff_ste(t,b) = .5*(quantile(off_rate_all,.75)-quantile(off_rate_all,.25));
            init_rate_all =  (pi2_vec.*initiation_rates(2,q_filter)'+...
                pi3_vec.*initiation_rates(3,q_filter)')./(pi2_vec+pi3_vec);  
            init_rate_eff_med(t,b) = median(init_rate_all);
            init_rate_eff_ste(t,b) = .5*(quantile(init_rate_all,.75)-quantile(init_rate_all,.25));
        end
    end
end
%%% Make Output Struct With Relevant Fields
hmm_results = struct;
for i = 1:length(bin_range)
    for t = 1:length(time_index)
        ind = (i-1)*length(time_index) + t;
%         hmm_results(ind).N = bin_counts(i);    
        hmm_results(ind).n_boots_used = n_boots_used(t,i);
        hmm_results(ind).n_boots_total = n_boots_total(t,i);
        hmm_results(ind).initiation_mean = med_initiation(:,i,t);
        hmm_results(ind).initiation_std = std_initiation(:,i,t);    
        hmm_results(ind).occupancy_mean = med_occupancy(:,i,t);
        hmm_results(ind).occupancy_std = std_occupancy(:,i,t);        
        hmm_results(ind).pi0_mean = avg_pi0(:,i,t);
        hmm_results(ind).pi0_std = std_pi0(:,i,t);
        hmm_results(ind).dwell_mean = med_dwell(:,i,t);
        hmm_results(ind).dwell_std = std_dwell(:,i,t);    
        hmm_results(ind).A_mean = avg_A(t,:,i);        
        hmm_results(ind).R_orig_mean = med_R_orig(t,:,i);
        hmm_results(ind).R_orig_std = std_R_orig(t,:,i);
        hmm_results(ind).R_fit_mean = med_R_fit(t,:,i);
        hmm_results(ind).R_fit_std = std_R_fit(t,:,i);
        hmm_results(ind).noise_mean = med_noise(i,t);
        hmm_results(ind).noise_std = std_noise(i,t);
        hmm_results(ind).N_dp = med_n_dp(i,t);
        hmm_results(ind).N_tr = med_n_tr(i,t);
%         hmm_results(ind).t_inf_effective = med_eff_time(i,t);
        hmm_results(ind).noise_mean = med_noise(i,t);
        hmm_results(ind).time_vec = time_vec;
        hmm_results(ind).bin_vec = bin_vec;
        hmm_results(ind).initiation_array = initiation_rates;
        hmm_results(ind).R_orig_array = R_orig_array;
        hmm_results(ind).R_fit_array = R_fit_array;
        if K == 3
            hmm_results(ind).on_ratio_mean = med_on_ratio(t,i);
            hmm_results(ind).on_ratio_ste = std_on_ratio(t,i);
            hmm_results(ind).off_ratio_mean = med_off_ratio(t,i);
            hmm_results(ind).off_ratio_ste = std_off_ratio(t,i);
            hmm_results(ind).init_ratio_mean = med_init_ratio(t,i);
            hmm_results(ind).init_ratio_ste = std_init_ratio(t,i);
            hmm_results(ind).eff_on_rate = on_rate_eff_med(t,i);
            hmm_results(ind).eff_on_ste = on_rate_eff_ste(t,i);
            hmm_results(ind).eff_off_rate = off_rate_eff_med(t,i);
            hmm_results(ind).eff_off_ste = off_rate_eff_ste(t,i);
            hmm_results(ind).eff_init_rate = init_rate_eff_med(t,i);
            hmm_results(ind).eff_init_ste = init_rate_eff_ste(t,i);
        end
        hmm_results(ind).binID = bin_range(i);
        hmm_results(ind).t_inf = time_index(t);
        hmm_results(ind).alpha = alpha;
        hmm_results(ind).dT = Tres;
        hmm_results(ind).fluo_type = fluo_type; 
        hmm_results(ind).clipped = clipped; 
        hmm_results(ind).clipped = clipped_ends; 
    end
end
save([OutPath '/hmm_results_t_window' num2str(t_window)  '_t_inf' num2str(t_inf)  '.mat'],'hmm_results')
disp(['compiled results saved to: ' OutPath])
