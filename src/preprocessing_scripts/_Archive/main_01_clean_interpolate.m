addpath('../utilities')
%---------------------Set Cleaning and Summary Parameters-----------------%
%Which Nuclear Cycles do we wish to include?
nuclear_cycles = [14];
%Incoming Time Resolution
Tres_in = 10.2;
%memory (best guess at memory)
w_in_appx = 14;
%Min dp per trace
min_dp_in = w_in_appx*2;
%outgoing mem
w_out_appx = 7;
%Min dp per trace
min_dp_out = w_out_appx*2;
%Interp Res
Tres_interp = 20;
%Length of MS2 loops in time steps
alpha_rat = 1302/6444;
%Designate sets to remove from final, cleaned set (due to QC issues)
rm_names = {'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\2017-05-14-A150umP_eve_5uW\CompiledParticles_classifier_2017_05_14_A150umP_eve.mat',...
            'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\2017-05-13-A200umP_eve_5uW\CompiledParticles_classifier_2017-05-14-A150umP_eve_5uW.mat',...
            'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\2017-07-06-A200umP_eve_5uW\CompiledParticles.mat',...
            'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\2017-07-10-P100umA_eve_5uW\CompiledParticles.mat',...
            'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\2017-07-09-A150umP_eve_5uW\CompiledParticles.mat'
            };
start_time = 0;
stop_time = 60;
%------------------------Define Project Vars------------------------------%
datatype = 'weka';
project = ['eve7stripes_inf_2017_09_25'];
%Output PNG files for each trace?
print_traces = 1;
%---------------------------Set Paths-------------------------------------%
inpath = ['../../dat/' project '/' 'raw_traces_' project '.mat'];
outpath = ['../../dat/' project '/'];
figpath = ['../../fig/' project '/'];

dataname = ['inference_traces_t' num2str(Tres_interp) '_' project '.mat'];
tracepath = [figpath '/traces_dT' num2str(Tres_interp) '/'];

if exist(figpath) ~= 7
    mkdir(figpath);
end
if exist(tracepath) ~= 7
    mkdir(tracepath);
end

%----------------Load Traces and Perform First Pass QC--------------------%
%Load raw traces (saved in struct titled "trace_struct")
load(inpath);
%remove flagged sets
rm_list = [];
for i = 1:length(trace_struct)
    for j = 1:length(rm_names)
        if strcmp(trace_struct(i).set,rm_names{j})
            rm_list = [rm_list i];
        end
    end
end
index_vec = 1:length(trace_struct);
trace_struct = trace_struct(~ismember(index_vec,rm_list));

%Determine time resolution of raw data
set_list = unique([trace_struct.setID]);
d_times = [];
trace_dt = [];
set_ids = [];
d_fig = figure;
hold on
for i = 1:length(set_list)
    subplot(4,4,i);    
    set_struct = trace_struct([trace_struct.setID]==set_list(i));
    dts = [];
    for t = 1:length(set_struct)
        t = trace_struct(i).time;
        dts = [dts diff(t)];
        d_times = [d_times diff(t)];
    end
    histogram(dts(dts<40),100)
    title(['Set ' num2str(set_list(i))])
end

%%
%Data quality checks. Discard suspect data points.
%Currently adjusting for the following issues:
% -"Sudden" drops to zeros or rises from zero. It is often the case that the first 1-4 steps
%       will be missed before fluo ramps up past a (variable) detection threshold
% -"Blips" where fluo drops and rises suddenly or vice versa. Only filtering for two step
%       cases. Unclear to me how best to address longer lasting anamolies automatically 
% -False "heads" and "tails"

%Perform cleaning seperately for each stripe
stripes = [trace_struct.stripe_id];
stripe_set = unique(stripes(~isnan(stripes)));
%set to pool all cleaned traces
trace_struct_final = [];
for s = stripe_set
    stripe_struct = trace_struct([trace_struct.stripe_id]==s);
    index_vec = 1:length(stripe_struct);
    %set scale for cleaning
    big_jump = prctile([stripe_struct.fluo],99); 
    %Remove low quality traces
    field_names = fieldnames(stripe_struct);
    %Structure to store traces after first round of cleaning
    stripe_struct_clean1 = [];
    flags = 0;
    %Remove infeasibly large rises and falls first. Save as much of the
    %remaining trace as possible
    for i = 1:length(stripe_struct)
        %Load full trace, including intervening NaN's
        trace = stripe_struct(i).fluo_full;
        
        %Exclude single points
        if length(trace) < 2        
            continue
        end    
        time = stripe_struct(i).time_full;
        trace_ind = 1:length(trace);
        %Null assumption is that all clusters of 3 or more NaNs are 0s. Single
        %or double NaNs are assumed to have been missed nonzero dps, will be interpolated
        trace_nans = isnan(trace);   
        %Look for clusters of 3 or more NaNs
        span = floor(60/Tres_in);
        nf = 0;
        cn = 0;
        nan_stretches = zeros(1,length(trace_nans));
        for t = 1:length(trace_nans)
            if trace_nans(t) == 1 && nf == 0
                nf = 1;
                cn = cn + 1;                
            elseif trace_nans(t) == 1 && nf == 1
                cn = cn + 1;               
            elseif trace_nans(t) == 0
                if nf == 1 && cn > span
                    nan_stretches(t-cn:t-1) = 1;
                end
                cn = 0;
                nf = 0;
            end
        end
                
        %Set clusters to zeros
        trace(nan_stretches==1) = 0;
        %interpolate remaining NaNs    
        query_times = time(isnan(trace));
        interp_t = time(~isnan(trace));
        interp_f = trace(~isnan(trace));
        new_f = interp1(interp_t,interp_f,query_times);
        trace(isnan(trace)) = new_f; 

        %deal with negative values
        trace(trace<0) = 0;
        %find zeros
        trace_zeros = trace==0;
        trace_nz = find(trace);
        %Look for stretches of 2+ zeros to land if a trace is broken
        trace_zeros = find(trace_zeros.*(0==[diff(trace_zeros) 0]));
        %Differentiate
        tr_d = diff(trace);    
        tr_dd = [0 diff(diff(trace)) 0]; 
        %Catch quick rises first 
        rise_points = find(([tr_d 0] > 8*big_jump/w_in_appx).*(trace==0).*(trace_ind~=1));
        problems = [];
        problems = [problems rise_points];
        %same story for drops but looking ahead
        fall_points = find(([0 tr_d] < -8*big_jump/w_in_appx).*(trace==0));    
        problems = [problems rise_points];
        problems = sort(problems);
        problems = [problems length(trace)];
        flags = flags + length(problems) -1;    
        trace_zeros = trace_zeros(~ismember(trace_zeros,problems));
        %Try to find contiguous stretches without anamalous rises and falls
        temp_struct = struct;
        start = 1;
        j_pass = 1;
        for j = 1:length(problems)
            trf = trace(start:problems(j));
            tif = time(start:problems(j));
            if length(trf) > 1
                temp_struct(j_pass).fluo_adjusted = trf;
                temp_struct(j_pass).time_adjusted = tif;
                j_pass = j_pass + 1;            
            end
            %If we find another start point, make sure it is far removed from
            %problem area
            s_options = trace_zeros(trace_zeros > problems(j));
            if isempty(s_options)
                break
            end
            %Find new start
            start = s_options(1);
        end
        if j_pass == 1
            continue
        end
        for j = 1:length(temp_struct)
            for f = 1:length(field_names)           
                temp_struct(j).(field_names{f}) = stripe_struct(i).(field_names{f});
            end
        end
        stripe_struct_clean1 = [stripe_struct_clean1 temp_struct];   
    end
    
    % Remove intermediate-sized Jumps. These will be smoothed via interpolation    
    adjust_flags = 0;
    st_flags = 0;
    stripe_struct_clean2 = [];
    for i = 1:length(stripe_struct_clean1)
        temp_struct = stripe_struct_clean1(i);    
        trace = stripe_struct_clean1(i).fluo_adjusted;
%         ap = stripe_struct_clean1(i).ap_vector;
%         xp = stripe_struct_clean1(i).xPos;
%         yp = stripe_struct_clean1(i).yPos;
        time = stripe_struct_clean1(i).time_adjusted;
        tr_d = diff(trace);    
        tr_dd = abs([0 diff(diff(trace)) 0]);
        trace_ind = 1:length(trace);
        rise_points = ([tr_d 0] > 4*big_jump/w_in_appx).*(trace==0).*(trace_ind~=1);
        fall_points = ([0 tr_d] < -4*big_jump/w_in_appx).*(trace==0);
        start_flag =  trace(1) > big_jump/w_in_appx;    
        %any large transient spike surrounded by nonzero values is being interpolated
        spike_points = ([0 abs(tr_d) > big_jump/w_in_appx]).*([abs(tr_d) 0] > big_jump/w_in_appx).*(tr_dd > 4*2*big_jump/w_in_appx);
        adjust_points = rise_points + fall_points + spike_points;
                                
        % if trace starts high, join leading zeros and interpolate
        if start_flag == 1
            trace = [0 0 trace];
            time = [time(1)-2*Tres_interp time(1)-Tres_interp time];
            adjust_points = [0 1 adjust_points];
        end    

        %Remove trouble points. will be interpolated
        trace = trace(adjust_points==0);
        time = time(adjust_points==0);        

%         temp_struct.time_adjusted_fluo = time_f;
        temp_struct.time_adjusted = time;
        temp_struct.fluo_adjusted = trace;
        
        stripe_struct_clean2 = [stripe_struct_clean2 temp_struct];
        adjust_flags = adjust_flags + sum(adjust_points==1);
    end
    index_vector = 1:length(stripe_struct_clean2);
    rm_list = [];

    for i = 1:length(stripe_struct_clean2)
        ff = stripe_struct_clean2(i).fluo;
        if length(ff) <= min_dp_in
            rm_list = [rm_list i];
        end
    end
    stripe_struct_final = stripe_struct_clean2(~ismember(index_vector,rm_list));
    trace_struct_final = [trace_struct_final stripe_struct_final];
end
%%
for i = 1:10:length(trace_struct_final)
    hold off
    plot(trace_struct_final(i).time_adjusted,trace_struct_final(i).fluo_adjusted,'-o','LineWidth',2)
    hold on
    plot(trace_struct_final(i).time,trace_struct_final(i).fluo,'-s')
    pause(1.5)
end
%% Interpolate to Achieve Desired Memory
interp_struct_raw = struct;
i_pass = 1;
for i = 1:length(trace_struct_final)
    fluo = trace_struct_final(i).fluo_adjusted;
    time = trace_struct_final(i).time_adjusted;
%     time_a = trace_struct_final(i).time_adjusted;
    
    fluo_orig = trace_struct_final(i).fluo;
    time_orig = trace_struct_final(i).time;    
    
    ap_orig = trace_struct_final(i).ap_vector;
    xPos_orig = trace_struct_final(i).xPos;
    yPos_orig = trace_struct_final(i).yPos;
    
    time_interp = min(time):Tres_interp:(min(time) + Tres_interp*floor((max(time)-min(time))/Tres_interp));
    fluo_interp = interp1(time,fluo,time_interp);
    if min(time_orig) > min(time_interp)        
        ap_orig = [ap_orig(1) ap_orig];
        xPos_orig = [xPos_orig(1) xPos_orig];
        yPos_orig = [yPos_orig(1) yPos_orig];
        time_orig = [time_interp(1) time_orig];
        fluo_orig = [NaN fluo_orig];
    end
    ap_interp = interp1(time_orig,ap_orig,time_interp);
    xPos_interp = interp1(time_orig,xPos_orig,time_interp);
    yPos_interp = interp1(time_orig,yPos_orig,time_interp);   
    if sum(isnan(ap_interp)) > 0
        error('Unbalanced Time Frames')
    end
    f_start = find(fluo_interp,1);
    f_stop = find(fluo_interp,1,'last');
    fluo_interp = fluo_interp(max(1,f_start-1):min(length(fluo_interp),f_stop+1));
    time_interp = time_interp(max(1,f_start-1):min(length(time_interp),f_stop+1));
    st_point = stop_time*60;
    if length(fluo_interp) >= min_dp_in
        interp_struct_raw(i_pass).fluo = fluo_interp(1:min(length(fluo_interp),st_point));
        interp_struct_raw(i_pass).time = time_interp(1:min(length(fluo_interp),st_point));
        interp_struct_raw(i_pass).fluo_orig = fluo_orig;
        interp_struct_raw(i_pass).time_orig = time_orig;
        
        interp_struct_raw(i_pass).ap_vector = ap_interp;
        interp_struct_raw(i_pass).xPos = xPos_interp;
        interp_struct_raw(i_pass).yPos= yPos_interp;
        field_names = fieldnames(trace_struct_final);
        for f = 1:length(field_names)
            fn = field_names{f};
            if ~strcmp(fn,'fluo')&&~strcmp(fn,'time')&&~strcmp(fn,'fluo_orig')&&~strcmp(fn,'time_orig')...
                && ~strcmp(fn,'xPos')&&~strcmp(fn,'yPos')&&~strcmp(fn,'ap_vector')
                interp_struct_raw(i_pass).(fn) = trace_struct_final(i).(fn);
            end
        end           
        interp_struct_raw(i_pass).N = length(fluo_interp);
        interp_struct_raw(i_pass).dT = Tres_interp;        
        interp_struct_raw(i_pass).w_cleaning = w_out_appx;
        interp_struct_raw(i_pass).alpha_frac = alpha_rat;
        interp_struct_raw(i_pass).start_time = start_time;
        interp_struct_raw(i_pass).stop_time = stop_time;
        i_pass = i_pass + 1;
    end
end

% define decision point = 50th prctile of nonzero pojnts
f_int = [interp_struct_raw.fluo];
f_int = f_int(f_int>0);
midpoint = prctile(f_int,50);
% check for infeasible spikes from zero level
problems = 0;
rejects = [];
problem_segments = [];
interp_struct = [];
for i = 1:length(interp_struct_raw)
    end_points = [];
    start_points = [];
    fluo = interp_struct_raw(i).fluo;
    d_f = abs(diff(fluo));
    zs = find(fluo==0);
    s = 1;
    for j = 1:(length(zs)-1)        
        %Find blips > 1 min and < twice memory (min ris/fall time)
        if zs(j+1) - zs(j) > floor(60/Tres_interp) && zs(j+1) - zs(j) < 2*w_out_appx
            % if a has a mean fluorescence greater than 50th percentile,
            % include in list of problem sections
            if mean(d_f(zs(j):zs(j+1)-1)) > midpoint
                problems = problems + 1;
                problem_segments = [problem_segments {fluo(zs(j):zs(j+1))}];
                %set new end point befor issue and start point right after
                end_points = [end_points zs(j)]; 
                start_points = [start_points zs(j+1)];
            end
        end
    end              
    start_points = [1 start_points];
    end_points = [end_points length(fluo)];
    % see what remains
    for ind = 1:length(start_points)  
        seg = fluo(start_points(ind):end_points(ind));
        if length(seg) >= min_dp_in
            temp = interp_struct_raw(i);
            temp.time = temp.time(start_points(ind):end_points(ind));
            temp.fluo = temp.fluo(start_points(ind):end_points(ind));
            temp.xPos = temp.xPos(start_points(ind):end_points(ind));
            temp.yPos = temp.yPos(start_points(ind):end_points(ind));
            temp.ap_vector = temp.ap_vector(start_points(ind):end_points(ind));
            if length(temp.time) ~= length(temp.ap_vector)
                error('wtf')
            end            
            interp_struct = [interp_struct temp];            
        else
            rejects = [rejects {seg}];
        end
    end
end
%Look  start and end zeros tails
rm_list = [];
for i = 1:length(interp_struct)
    trace = interp_struct(i).fluo;   
    rm_points = zeros(size(trace));
    if max(trace(1:3)) ~= 0 && sum(trace(4:8)) == 0
        trace(1:3) = 0;
        new_start = find(trace,1);        
        rm_points(1:new_start-1) = 1; 
    end
    if max(trace(end-2:end)) ~= 0 && sum(trace(end-7:end-3)) == 0
        trace(end-2:end) = 0;
        new_end = find(trace,1,'last');        
        rm_points(new_end+1:end) = 1; 
    end
    new_fluo = trace(rm_points==0);
    new_time = interp_struct(i).time(rm_points==0);
    new_ap = interp_struct(i).ap_vector(rm_points==0);
    new_x = interp_struct(i).xPos(rm_points==0);
    new_y = interp_struct(i).yPos(rm_points==0);
    % One last filter to find first and last non-zero points
    start = max(1,find(new_fluo,1)-1);
    stop = min(length(new_fluo),find(new_fluo,1,'last')+1);
    interp_struct(i).fluo = new_fluo(start:stop);
    interp_struct(i).time = new_time(start:stop);
    interp_struct(i).ap_vector = new_ap(start:stop);
    interp_struct(i).xPos = new_x(start:stop);
    interp_struct(i).yPos = new_y(start:stop);
    interp_struct(i).OrigParticleLong = repelem(interp_struct(i).OriginalParticle,length(new_y(start:stop)));
    interp_struct(i).StripeIDLong = repelem(interp_struct(i).stripe_id,length(new_y(start:stop)));
    interp_struct(i).StripeSubIDLong = repelem(interp_struct(i).stripe_sub_id,length(new_y(start:stop)));
    interp_struct(i).SetIDLong = repelem(interp_struct(i).setID,length(new_y(start:stop)));
    interp_struct(i).StripeCenterLong = repelem(interp_struct(i).stripe_center_ap,length(new_y(start:stop)));
    
    if length(interp_struct) < min_dp_in || max(new_fluo) == 0
        rm_list = [rm_list i];
    end
end
ind_vec = 1:length(interp_struct);
interp_struct = interp_struct(~ismember(ind_vec,rm_list));

% Save data
save([outpath dataname],'interp_struct');
% Generate csv
output = [[interp_struct.SetIDLong]',[interp_struct.OrigParticleLong]', ...
         [interp_struct.StripeIDLong]', [interp_struct.StripeSubIDLong]',...
          [interp_struct.time]', [interp_struct.fluo]',...
          [interp_struct.StripeCenterLong]', [interp_struct.ap_vector]', ...
          [interp_struct.xPos]', [interp_struct.yPos]'];
    
header = {'SetID','ParticleID','StripeID','StripeSubID','Seconds', 'Fluo',...
    'StripeCenter','AP','x','y'};
csvwrite_with_headers([outpath '\' dataname '_dT' num2str(Tres_interp) '_longform.csv'], ...
                       output, header); 
%% Save Cleaned Set of un-interpolated traces
%Look  start and end zeros tails
rm_list = [];
for i = 1:length(trace_struct)
    n = length(trace_struct(i).fluo);
    % One last filter to find first and last non-zero points
    trace_struct(i).OrigParticleLong = repelem(trace_struct(i).OriginalParticle,n);
    trace_struct(i).StripeIDLong = repelem(trace_struct(i).stripe_id,n);
    trace_struct(i).StripeSubIDLong = repelem(trace_struct(i).stripe_sub_id,n);
    trace_struct(i).SetIDLong = repelem(trace_struct(i).setID,n);
    trace_struct(i).StripeCenterLong = repelem(trace_struct(i).stripe_center_ap,n);
    
    if length(trace_struct) < min_dp_in || max(trace_struct(i).fluo) == 0
        rm_list = [rm_list i];
    end
end
ind_vec = 1:length(trace_struct);
trace_struct = trace_struct(~ismember(ind_vec,rm_list));


% Generate csv
output = [[trace_struct.SetIDLong]',[trace_struct.OrigParticleLong]', ...
         [trace_struct.StripeIDLong]', [trace_struct.StripeSubIDLong]',...
          [trace_struct.time]', [trace_struct.fluo]',...
          [trace_struct.StripeCenterLong]', [trace_struct.ap_vector]', ...
          [trace_struct.xPos]', [trace_struct.yPos]'];
    
header = {'SetID','ParticleID','StripeID','StripeSubID','Seconds', 'Fluo',...
    'StripeCenter','AP','x','y'};
csvwrite_with_headers([outpath '\' project '_raw_longform.csv'], ...
                       output, header); 
%--------------------Print Sample Traces----------------------------------%
%% If desired, save individual trace plots
if print_traces
    for i = 1:10:length(interp_struct)
        t_fig = figure('Visible','off');
        hold on
        plot(interp_struct(i).time_orig / 60, interp_struct(i).fluo_orig, '-o','Linewidth',1.5)
        plot(interp_struct(i).time / 60, interp_struct(i).fluo,'-s', 'Linewidth',1.5)
        title(['Original vs. Final: Trace ' num2str(i)])
        legend('Raw', 'Interpolated');
        xlabel('Minutes');
        grid on        
        saveas(t_fig, [tracepath, 'trace_' num2str(i) '.png'],'png');
    end
end
