%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = '../../../data/';
outpath = '../../processed_data/';
outName = 'eveSet_2017_06_21_250_only.mat';
if exist(outpath) ~= 7
    mkdir(outpath);
end
files = dir(folder_path);

filenames = {};
for i = 1:length(files)
    if ~isempty(strfind(files(i).name,'.mat')) && ~isempty(strfind(files(i).name,'250'))
        filenames = [filenames {files(i).name}];
    end
end

traces_no_nan = struct;
traces_interp = struct;
diff_traces = 0;
%Data structure to store extracted trace sets

trace_struct = struct;
i_iter = 1;
for k = 1:length(filenames)
    raw_data  = load([folder_path filenames{k}]);
    time = raw_data.ElapsedTime*60;
    time = time - time(raw_data.nc14);
    traces = raw_data.AllTracesVector;
    nc14_filter = raw_data.ncFilter(:,end);
    traces_14 = traces(:,nc14_filter == 1);
    %Set Nans to -1
    traces_14(traces_14==0) = -1;
    traces_14(isnan(traces_14)) = 0;
    for i = 1:size(traces_14,2)
        raw_trace = traces_14(:,i);
        start = find(raw_trace,1);
        stop = find(raw_trace,1,'last');
        trunc_trace = [raw_trace(start:stop)'];
        trunc_time = time(start:stop); 
        [~, apPos] = max(raw_data.APFilter(i,:));
        
        trace_struct(i_iter).fluo = trunc_trace;
        trace_struct(i_iter).time = trunc_time;
        trace_struct(i_iter).AP = apPos;
        trace_struct(i_iter).set = filenames{k};
        i_iter = i_iter + 1;
    end
end

% Run data quality checks. Discard suspect dps. Will be replaced in interp step

%Get 90th percentile for point-to-point deltas
ref_len = prctile(abs(diff([trace_struct.fluo])),90);
adjustments = 0;
for i = 1:length(trace_struct)
    trace = trace_struct(i).fluo;
    tr_d = [0 diff(trace)];
    tr_dd = [0 diff(diff(trace)) 0];
    rm_list = [];
    for j = 1:length(trace)
        if abs(tr_d(j)) > .5 * ref_len && trace(j) == 0
            rm_list = [rm_list j];
            adjustments = adjustments +1;
        elseif abs(tr_dd) > 1.5*ref_len
            rm_list = [rm_list j];
            adjustments = adjustments +1;
        end
    end
    trace_struct(i).fluo = trace(~ismember(1:length(trace), rm_list));
    trace_struct(i).time = trace_struct(i).time(~ismember(1:length(trace),rm_list));
end
        
% Interpolate to Achieve Tractable Memory
%Define Desired res and associated param
mem = 8;
T_elong = 150;
Tres = T_elong / mem;
%Set minimum trace length (in time steps)
min_len = mem;
i_pass = 1;
interp_struct = struct;
for i = 1:length(trace_struct)
    fluo = trace_struct(i).fluo;
    time = trace_struct(i).time;
    
    time_interp= min(time) + Tres*(0:floor((max(time)-min(time))/Tres));
    fluo_interp = interp1(time,fluo,time_interp);
    if length(time_interp) > min_len
        interp_struct(i_pass).fluo = fluo_interp;
        interp_struct(i_pass).time = time_interp;
        interp_struct(i_pass).fluo_orig = fluo;
        interp_struct(i_pass).time_orig = time;
        interp_struct(i_pass).AP = trace_struct(i).AP;
        interp_struct(i_pass).set = trace_struct(i).set;
        interp_struct(i_pass).N = length(fluo_interp);
        interp_struct(i_pass).dT = Tres;
        interp_struct(i_pass).T_elong = T_elong;
        interp_struct(i_pass).w = mem;
        interp_struct(i_pass).alpha = (60/204)*mem;
        i_pass = i_pass + 1;
    end
end

% Truncate Traces with Unreasonably Long Pauses
%Set reference p_off
p_off_ref = .9;
p_cp = 1;
power = 0;
while p_cp > .005;
    power = power + 1;
    p_cp = p_cp * p_off_ref;
end

cuts = 0;
rm_ids = [];
for i = 1:length(interp_struct)
    fluo = interp_struct(i).fluo;
    time = interp_struct(i).time;
    FluoBin = fluo==0;
    ct = 0;
    for j = 1:length(fluo)
        ct = ct*FluoBin(j) + FluoBin(j);
        if ct > power
            fluo_cut = fluo(1:j-ct);
            time_cut = time(1:j-ct);
            cuts = cuts + 1;
            if length(fluo_cut) >= min_len
                interp_struct(i).fluo = fluo_cut;
                interp_struct(i).time = time_cut;
                interp_struct(i).N = length(fluo_cut);
            else
                rm_ids = [rm_ids i];
            end
            break
        end
    end
end
if ~isempty(rm_ids)
    interp_struct = interp_struct(~ismember(1:length(interp_struct),rm_ids));
end
display([num2str(cuts) ' traces truncated due to infeasibly long pauses']);

save([outpath outName], 'interp_struct');