function gillespie = synthetic_rate_gillespie(seq_length, alpha, ...
                                K, w, R, deltaT, r_emission, noise, pi0)

    % Generates a fluorescence sequence for the given model parameters
    % using the Gillespie algorithm.
    % 
    % INPUTS
    % seq_length: length of the sequence in time steps
    % alpha: length of the MS2 loop in time steps
    % K: number of naive states
    % w: memory
    % R: transition rate matrix [sec^{-1}]
    % deltaT: time resolution [sec]
    % r_emission: emission rates [sec^{-1}]
    % noise: Gaussian noise [a.u.]
    % pi0: initial state probability mass function (pmf)
    % 
    % OUTPUTS
    % gillespie: structure that contains the synthetic data info
    %   .fluo: synthetic trace generated without taking the MS2 loops
    %          into account
    %   .fluo_MS2: synthetic trace generated by taking the MS2 loops 
    %              into account
    %   .fluo_no_noise: synthetic trace generated without taking the MS2 
    %                   loops into account and with no Gaussian noise
    %   .fluo_MS2_no_noise: synthetic trace generating by taking the MS2 
    %                       loops into account but with no Gaussian noise
    %   .transition_times: times when the transitions happened
    %   .transition_count: number of transitions in the simulated trace
    %   .naive_states: sequence of the naive states of the system at time
    %                  moments specified in the variable 'times'
    
    % uniformly distributed time points with resolution deltaT 
    % and length seq_length
    times_unif = (1:seq_length) * deltaT;
    
    % fluorescence sequence at time points 'time_unif' when the effect of
    % the MS2 loops is taken into account
    fluo_MS2 = zeros(1, seq_length);
    
    % arrays for storing the naive states, times of transitions and 
    % cumulative fluorescence values during the Gillespie simulation
    naive_states = [];
    times = [];
    fluo_cum = [];
    
    % duration of the simulated process
    t_max = seq_length * deltaT;
    
    % first state obtained by sampling the initial state pmf
    naive_states(1) = randsample(1:K, 1, true, pi0);
    
    % assign 0 values to the time and fluo_cum arrays at t = 0
    times(1) = 0;
    fluo_cum(1) = 0;
    
    % variable to keep track of the current reaction time
    t = 0;
    
    % indexing variable to keep track of the current transition
    ind = 1;
    
    % generate a sequence of naive states using the Gillespie algorithm
    while (t < t_max)
        
        lambda = -R(naive_states(ind),naive_states(ind));
        dt = exprnd(1/lambda);
        t = t + dt;
        
        rates = R(:,naive_states(ind));
        rates(naive_states(ind)) = 0;
        
        probs = rates / lambda;
        
        naive_states(ind + 1) = randsample(1:K, 1, true, probs);
        ind = ind + 1;
        times(ind) = t;
        
        fluo_cum(ind) = fluo_cum(ind-1) + (times(ind) - times(ind-1)) * ...
            r_emission(naive_states(ind - 1));
    end
    
    % number of transitions in the simulated trace
    transition_count = ind-1;
    
    % find the fluorescence at the uniformly distributed points without
    % accounting for the MS2 loops
    fluo_cum_interp = interp1(times, fluo_cum, times_unif);
    fluo_cum_interp_shift = [zeros([1, w]), fluo_cum_interp(1:(end-w))];
    fluo = fluo_cum_interp - fluo_cum_interp_shift;
    
    % find the fluorescence at the uniformly distributed points taking the
    % MS2 loops into account
    for k = 1:seq_length
        t_end = times_unif(k);
        t_start = max([0, t_end - w*deltaT]);
        
        ind_start_before = find(times <= t_start);
        i_start = ind_start_before(end);
        
        ind_end_before = find(times <= t_end);
        i_end = ind_end_before(end);
        
        times_window = times(i_start:(i_end+1));
        times_window(1) = t_start;
        times_window(end) = t_end;
        
        naive_window = naive_states(i_start:(i_end+1));
        
        for i = 1:(length(naive_window)-1)
            
            t1 = t_end - times_window(i+1);
            t2 = t_end - times_window(i);
            fluo_MS2(k) = fluo_MS2(k) + r_emission(naive_window(i)) * ...
                          ms2_loading_coeff_integral(alpha, w, deltaT, t1, t2);
            
            % old approach
%             fluo_MS2(k) = fluo_MS2(k) + r_emission(naive_window(i)) * ...
%                 ms2_loading_coeff_integral(alpha, w, deltaT, ...
%                     times_window(i)-t_start, times_window(i+1)-t_start);
        end
    end

    % Gaussian noise
    gauss_noise = normrnd(0,noise,1,seq_length);
    
    % add a Gaussian noise
    fluo_noise = fluo + gauss_noise;
    fluo_MS2_noise = fluo_MS2 + gauss_noise;
    
    % the output structure
    gillespie = struct('fluo', fluo_noise, 'fluo_MS2', fluo_MS2_noise, ...
        'fluo_no_noise', fluo, 'fluo_MS2_no_noise', fluo_MS2, ...
        'transition_times', times, 'transition_count', transition_count, ...
        'naive_states', naive_states);