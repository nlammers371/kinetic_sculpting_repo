%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\DropboxSingleTraces\';
sub_paths = {'Eve2_ML', 'Eve2_orig'};
key_names = {'ML', 'ORIG'};
project = 'mHMMeve2_ml_comparisons';
outpath = ['../projects/' project '/' ];
% Keyword to ensure only sets from current project are pulled
keyword = 'eve2_20sec_';
exclude = 'eve2_20sec_5';
if exist(outpath) ~= 7
    mkdir(outpath);
end
%%
meta_struct = struct;
for h = 1:2
    dir_struct = struct;
    i_pass = 1;
    dirinfo = dir([folder_path sub_paths{h}]);
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    subdirinfo = cell(length(dirinfo));
    for K = 1 : length(dirinfo)
      thisdir = dirinfo(K).name;
      % skip files lacking prject keyword
      if isempty(strfind(thisdir,keyword)) || ~isempty(strfind(thisdir,exclude))
          continue
      end
      subdir_struct = dir(fullfile(folder_path,sub_paths{h},thisdir, '*.mat'));
      if length(subdir_struct) > 5
          if i_pass == 1
              dir_struct(i_pass).files = subdir_struct;
              i_pass = 2;
          else
              dir_struct(i_pass).files = subdir_struct;
              i_pass = i_pass + 1;
          end
      end
    end

    filenames = {};
    for i = 1:length(dir_struct)
        subdir_struct = dir_struct(i).files;
        for j = 1:length(subdir_struct)
            if strfind(subdir_struct(j).name,'CompiledParticles') > 0
                filenames = [filenames {[subdir_struct(j).folder '\' subdir_struct(j).name]}];
            end
        end
    end

    traces_no_nan = struct;
    diff_traces = 0;
    %Data structure to store extracted trace sets
    trace_struct = struct;
    i_iter = 1;
    for k = 1:length(filenames)
        raw_data  = load([filenames{k}]);
        time = raw_data.ElapsedTime*60;
        time = time - time(raw_data.nc14);
        traces = raw_data.AllTracesVector;
        nc14_filter = raw_data.ncFilter(:,end);
        index_vector = 1:length(nc14_filter);
        ID14 = index_vector(nc14_filter == 1);
        %Set Nans to 0, true zeros to -1
        traces(traces==0) = -1000;
        traces(isnan(traces)) = 0;
        particles = raw_data.CompiledParticles;
        meta_struct(k).(['particles_' key_names{h}]) = particles;
        for i = ID14
            raw_trace = traces(:,i);
            nucleus = particles(i).Nucleus;
            particle = particles(i).OriginalParticle;
            % Skip single points
            if length(raw_trace(raw_trace ~= 0)) < 2
                continue
            end
            start = find(raw_trace,1);
            stop = find(raw_trace,1,'last');
            raw_trace(raw_trace==0) = nan;
            raw_trace(raw_trace==-1000) = 0;
            trunc_trace = [raw_trace(start:stop)'];
            trunc_time = time(start:stop); 
            [~, apPos] = max(raw_data.APFilter(i,:));
            trace_struct(i_iter).Index = particles(i).Index;
            trace_struct(i_iter).fluo = trunc_trace;
            trace_struct(i_iter).Frame = particles(i).Frame;
            trace_struct(i_iter).time = trunc_time;
            trace_struct(i_iter).AP = apPos;
            trace_struct(i_iter).set = filenames{k};
            trace_struct(i_iter).setID = k;
            trace_struct(i_iter).ncID = nucleus;
            trace_struct(i_iter).pID = particle;
            i_iter = i_iter + 1;
        end
    end
    % Run data quality checks. Discard suspect dps. Will be replaced in interp step
    % Get 95th percentile for point-to-point deltas
    diffs = abs(diff([trace_struct.fluo]));
    ref_len = prctile(diffs(~isnan(diffs)),99);
    adjustments = 0;
    for i = 1:length(trace_struct)
        trace = trace_struct(i).fluo;
        trace(isnan(trace)) = 0;
        tr_d = [0 diff(trace)];
        tr_dd = [0 diff(diff(trace)) 0];
        rm_list = [];
        for j = 1:length(trace)
            % remove "large" drops to zero
            if abs(tr_d(j)) > ref_len && trace(j) == 0
                rm_list = [rm_list j];
                adjustments = adjustments +1;
            % remove "large" transient spikes
            elseif abs(tr_dd) > 1.5*ref_len
                rm_list = [rm_list j];
                adjustments = adjustments +1;
            end
        end
        trace_struct(i).fluo = trace(~ismember(1:length(trace), rm_list));
        trace_struct(i).time = trace_struct(i).time(~ismember(1:length(trace),rm_list));
    end

    % Interpolate to Achieve Desired Memory
    % Define Desired res and associated param
    mem = 8;
    T_elong = 160;
    Tres = T_elong / mem;
    %Set minimum trace length (in time steps)
    min_len = mem;

    i_pass = 1;
    % define time grid. All traces will be shifted to conform to same grid
    t_start = Inf;
    for i = 1:length(trace_struct)
        if isempty(strfind(trace_struct(i).set,'eve2_20sec_3'))
            t_start  = min(t_start,floor(min(trace_struct(i).time)));
        end
    end
    t_stop = floor(max([trace_struct.time]));
    time_grid = t_start + Tres*(0:floor(t_stop/Tres));
    interp_struct = struct;
    for i = 1:length(trace_struct)
        if ~isempty(strfind(trace_struct(i).set,'eve2_20sec_3'))
            offset = t_start;
        else
            offset =0;
        end
        fluo = trace_struct(i).fluo;
        time = offset + trace_struct(i).time;

        time_interp = time_grid(time_grid >= min(time));
        time_interp = time_interp(time_interp <= max(time));
        fluo_interp = interp1(time,fluo,time_interp);
        if length(time_interp) > min_len
            interp_struct(i_pass).fluo = fluo_interp;
            interp_struct(i_pass).time = time_interp;
            interp_struct(i_pass).fluo_orig = fluo;
            interp_struct(i_pass).time_orig = time;
            interp_struct(i_pass).AP = trace_struct(i).AP;
            interp_struct(i_pass).set = trace_struct(i).set;
            interp_struct(i_pass).ncID = trace_struct(i).ncID;
            interp_struct(i_pass).pID = trace_struct(i).pID;
            interp_struct(i_pass).setID = trace_struct(i).setID;
            interp_struct(i_pass).Frame = trace_struct(i).Frame;
            interp_struct(i_pass).Index = trace_struct(i).Index;
            interp_struct(i_pass).N = length(fluo_interp);
            interp_struct(i_pass).dT = Tres;
            interp_struct(i_pass).T_elong = T_elong;
            interp_struct(i_pass).w = mem;
            interp_struct(i_pass).alpha = (60/204)*mem;
            i_pass = i_pass + 1;
        end
        meta_struct(k).(['interp_' key_names{h}]) = interp_struct;
    end
end
%%
orig_interp = meta_struct.interp_ORIG;
ml_interp = meta_struct.interp_ML;
%%
orig_nc = [orig_interp.ncID];
orig_set = [orig_interp.setID];
t_pass = 1;
for i = 1:length(ml_interp)
    nucleus = ml_interp(i).ncID;
    set = ml_interp(i).setID;
    if sum((orig_set == set).*(orig_nc == nucleus)) == 1
        orig_ind = find((orig_set == set).*(orig_nc == nucleus));
    	of = orig_interp(orig_ind).fluo_orig;
        ot = orig_interp(orig_ind).time_orig;
        mf = ml_interp(i).fluo_orig;
        mt = ml_interp(i).time_orig;
        mp = ml_interp(i).pID;
        op = orig_interp(i).pID;
        t_fig = figure('Visible','off','Position',[0 0 1024 512]);
        subplot(1,2,1)
        plot(ot, of,'-o', 'Linewidth',1.5)
        title(['Trace ' num2str(i) ' Set: ' num2str(set) ' Nuclues: ' num2str(nucleus) ...
            ' Particle: ' num2str(op) ' (Original)'])
        grid on
        axis([0 max([ot mt]) 0 1.1*max([of mf])]);
        subplot(1,2,2)
        plot(mt, mf,'-o', 'Linewidth',1.5)
        title(['Trace ' num2str(i) ' Set: ' num2str(set) ' Nuclues: ' num2str(nucleus) ...
            ' Particle: ' num2str(mp) ' (Weka)'])
        grid on
        axis([0 max([ot mt]) 0 1.1*max([of mf])]);
        saveas(t_fig, [outpath,'/traces/', 'trace_' num2str(i) '.png'],'png');
    end
end

