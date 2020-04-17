% Script to Compile Data Sets and Find Stripe Centers
% close all
% clear 
%------------------------Set Path Specs, ID Vars------------------------%
FolderPath = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\';
project = 'eve7stripes_inf_2018_03_27_final'; %Project Identifier
% folders
fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../../dat/' project '/']; % data mat directory
% fig subfolders
ap_pos_path = [fig_path 'ap_positioning/'];
fluo_path = [fig_path 'fluo_stats/'];
trace_name = [data_path 'raw_traces_' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_' project]; % names for compiled elipse struct

% make filepaths
mkdir(data_path);
mkdir([ap_pos_path '/stripe_fits']);
mkdir(fluo_path);
% cleaning params
keyword = '_30uW_550V'; % Keyword to ensure only sets from current project are pulled
xDim = 1024;
yDim = 256;
snippet_size = 15; % particles within snippet/2+1 are of concern
pre_post_padding = 10; % max mun frames for which nucleus can be MIA at start or end
% store set names
dirinfo = dir(FolderPath);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
cp_filenames = {}; % particles
ap_filenames = {}; % ap info
nc_filenames = {}; % nuclei
set_nums = [];
for d = 1 : length(dirinfo)
    thisdir = dirinfo(d).name;
    % Skip files lacking project keyword 
    if isempty(strfind(thisdir,keyword)) 
        continue
    end    
%     set_num_start_ind = strfind(thisdir,'_');
%     set_num_start_ind = set_num_start_ind(end);
%     set_num = str2num(thisdir(set_num_start_ind+1:end));        
%     set_nums = [set_nums set_num];
    % append file paths
    cp_name = dir([FolderPath '/' thisdir '/CompiledParticles*']);
    cp_name = cp_name(1).name;
    cp_filenames = [cp_filenames {[thisdir '/' cp_name]}];    
    ap_filenames = [ap_filenames {[thisdir '/APDetection.mat']}];    
    nc_filenames = [nc_filenames {[thisdir '/' thisdir '_lin.mat']}];           
end

trace_struct = struct; % Generate data structure to store extracted trace sets
schnitz_struct = []; % structure to store nucleis info

%%% compile traces and nuclei
j_pass = 0; 
total_matched = 0;
for i = 1:length(cp_filenames) % Loop through filenames    
    % read in raw files
    load([FolderPath ap_filenames{i}]) % AP Info   
    load([FolderPath nc_filenames{i}]) % Ellipse Info
    raw_data = load([FolderPath cp_filenames{i}]); % Particles    
    SetID = i;    
    % get angle between the x-axis and the AP-axis 
    APAngle = round(atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)))*360 / (2*pi));    
    %Correction for if APAngle is in quadrants II or III
    if coordPZoom(1)-coordAZoom(1) < 0       
        APAngle = APAngle + 180;
    end    
    %%% pull trace and nuclei variables
    time_raw = raw_data.ElapsedTime*60; % time vector            
    traces_raw = raw_data.AllTracesVector; % Array with a column for each trace    
    frames_raw = 1:length(time_raw); % Frame list    
    first_frame = raw_data.nc14; % Get frame that marks start of nc14
    last_frame = frames_raw(end);
    % Filter 
    traces_clean = traces_raw(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc14
    frames_clean = frames_raw(first_frame:end);     
    % Compile nucleus info
    s_cells = struct;
    e_pass = 1;    
    for e = 1:length(schnitzcells)
        e_frames = schnitzcells(e).frames;
        nuc_frames = e_frames(ismember(e_frames,frames_clean));
        if length(nuc_frames) > 2% skip nuclei not in nc14
            x = schnitzcells(e).cenx;
            x = x(ismember(nuc_frames,frames_clean));
            y = schnitzcells(e).ceny;
            y = y(ismember(nuc_frames,frames_clean));            
            s_cells(e_pass).xPos = x;
            s_cells(e_pass).yPos = y; 
            s_cells(e_pass).set_vector = repelem(SetID,length(x)); 
            s_cells(e_pass).APAngle = APAngle;
            % used to assess eligibility for on/off time analyses
            %Will be set to particle position for nuclei with matching
            %particle
            s_cells(e_pass).xPosParticle = NaN;
            s_cells(e_pass).yPosParticle = NaN;
            s_cells(e_pass).frames = nuc_frames';
            s_cells(e_pass).N = length(nuc_frames);
            s_cells(e_pass).Nucleus = e;                        
            s_cells(e_pass).ncID = eval([num2str(SetID) '.' iIndex(e,4)]);
            %Will be set to mean particle position for nuclei with matiching
            %particles
            s_cells(e_pass).xMean = mean(x);
            s_cells(e_pass).yMean = mean(y);
            s_cells(e_pass).time = time_clean(ismember(frames_clean,nuc_frames));                        
            s_cells(e_pass).setID = SetID;
            e_pass = e_pass + 1;
        end
    end
    %%% Now Particles
    e_index = [s_cells.Nucleus]; % Index vector to cross-ref w/ particles        
    fn = cp_filenames{i}; % Get filename to store in struct           
    % iterate through traces
    particles = raw_data.CompiledParticles;    
    j_init = j_pass;
    for j = 1:size(traces_clean,2)        
        raw_trace = traces_clean(:,j);     
        trace_start = find(~isnan(raw_trace),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs
        trace_full = raw_trace(trace_start:trace_stop)';        
        % skip small fragments      
        if length(raw_trace(~isnan(raw_trace))) < 3            
            continue
        end
%         elseif length(raw_trace(~isnan(raw_trace))) < 5            
%             short_flag = 1;
%         end
        j_pass = j_pass + 1;
        % Generate full length time and frame vectors     
        time_full = time_clean(trace_start:trace_stop);                
        frames_full = frames_clean(trace_start:trace_stop);
        % Pull variables from particle structure
        % Why the hell is there an APPos field with negative values?
        ap_positions_raw = particles(j).APpos;
        fov_xPos_raw = particles(j).xPos;
        fov_yPos_raw = particles(j).yPos;        
        pt_frames = particles(j).Frame;
        % Only take values from CP struct corresponding to frames in
        % filtered frames              
        ap_positions = ap_positions_raw(ismember(pt_frames,frames_full));        
        xPos = fov_xPos_raw(ismember(pt_frames,frames_full));        
        yPos = fov_yPos_raw(ismember(pt_frames,frames_full));        
        pt_frames = pt_frames(ismember(pt_frames,frames_full));
        % look for edge issues
        edge_frames = ((xDim-xPos) <= snippet_size)|(xPos <= snippet_size)|...
                          ((yDim-yPos) <= snippet_size)|(yPos <= snippet_size);
        trace_struct(j_pass).edge_flag = (sum(edge_frames)>0);        
        % Record info in trace struct        
        trace_struct(j_pass).cp_frames = pt_frames; % Keep track of frame correspondence
        trace_struct(j_pass).all_frames = frames_full;
        trace_struct(j_pass).nc14 = first_frame;
        trace_struct(j_pass).last_frame = last_frame;        
        trace_struct(j_pass).APAngle = APAngle;
        trace_struct(j_pass).xPos = xPos;
        trace_struct(j_pass).yPos = yPos;
        trace_struct(j_pass).ap_vector = ap_positions;
        trace_struct(j_pass).set_vector = repelem(SetID,length(ap_positions));
        trace_struct(j_pass).fluo = trace_full;        
        trace_struct(j_pass).time = time_full;
        trace_struct(j_pass).FluoError = particles(j).FluoError; %Estimated error in bkg subtraction     
        % For nuclei corresponding to particles, add fluo info and revise
        % mean position to align with particle pos
        nucleus = particles(j).Nucleus; 
        % Identifier variables        
        trace_struct(j_pass).Nucleus = nucleus;
        ncID = eval([num2str(SetID) '.' iIndex(nucleus,4)]);
        trace_struct(j_pass).ncID = ncID;
        particle = particles(j).OriginalParticle;        
        pID = eval([num2str(SetID) '.' iIndex(particle,4)]);                    
        trace_struct(j_pass).ParticleID = pID;
        trace_struct(j_pass).setID = SetID;
        trace_struct(j_pass).source_path = fn;                        
    end
    % Add Particle Info to Nuclei
    for j = j_init+1:j_pass
        Nucleus = trace_struct(j).Nucleus;
        nc_ind = find(e_index==Nucleus);
        if length(nc_ind) ~= 1
            error('Error: Problem with Particle-Nucleus Crossref')
        end
        total_matched = total_matched + 1;
        s_cells(nc_ind).xPosParticle = trace_struct(j).xPos;
        s_cells(nc_ind).yPosParticle = trace_struct(j).yPos;
        s_cells(nc_ind).xMean = mean(trace_struct(j).xPos);
        s_cells(nc_ind).yMean = mean(trace_struct(j).yPos);                
        s_cells(nc_ind).ParticleID = trace_struct(j).ParticleID;                
    end
    schnitz_struct = [schnitz_struct  s_cells];        
end
for i = 1:length(schnitz_struct)
    if isempty(schnitz_struct(i).ParticleID)
        schnitz_struct(i).ParticleID = NaN;
    end
end

%%% Look for trace fragments that belong together. Stitch them up

%%% This is premised on the trace start times being sorted in ascending
%%% order!
set_vec = [trace_struct.setID];
set_index = unique(set_vec);
% stitch together overlaps
remove_indices = [];
match_indices = [];
dupe_indices = [];
for s = 1:length(set_vec)
    set_struct_tr = trace_struct(set_index==set_vec(s));     
    base_indices = find(set_index==set_vec(s));
    start_t = [];
    stop_t = [];
    start_x = [];
    stop_x = [];
    start_y = [];
    stop_y = [];
    for i = 1:length(set_struct_tr)
        time = set_struct_tr(i).time;
        xPos = set_struct_tr(i).xPos;
        yPos = set_struct_tr(i).yPos;
        % add start and stop info
        start_t = [start_t time(1)];
        stop_t = [stop_t time(end)];
        start_x = [start_x xPos(1)];
        stop_x = [stop_x xPos(end)];
        start_y = [start_y yPos(1)];
        stop_y = [stop_y yPos(end)];
    end
    %%% enforce ascending sort order
    [start_t, si] = sort(start_t);
    stop_t = stop_t(si);
    start_x = start_x(si);
    stop_x = stop_x(si);
    start_y = start_y(si);
    stop_y = stop_y(si);
    base_indices = base_indices(si);
    if sum(diff(start_t)<0)>1
        error('Non-ascending start time order')
    end
    t_mat = repmat(start_t,length(start_t),1) - repmat(stop_t',1,length(stop_t));
    x_mat = repmat(start_x,length(start_x),1) - repmat(stop_x',1,length(stop_x));
    y_mat = repmat(start_y,length(start_y),1) - repmat(stop_y',1,length(stop_y));
    logic_mat = (sqrt((x_mat.^2 + y_mat.^2))<=5)&(t_mat>0)&(t_mat<120);
    logic_mat(eye(length(stop_y))==1) = 0; % remove diagonals
    overlap_ids = find(logic_mat) ;    
    tr_col = floor((overlap_ids-1)/length(stop_x))+1;
    tr_row = overlap_ids - (tr_col-1)*length(stop_x);
    tr_col = base_indices(tr_col);
    tr_row = base_indices(tr_row);        
    if length(unique(tr_col))~=length(tr_col) || length(unique(tr_row))~=length(tr_row)           
        warning('Duplicate Trace Fragments Detected. Removing.')
        col_mat = repmat(tr_col,length(tr_col),1)==repmat(tr_col',1,length(tr_col));
        row_mat = repmat(tr_row,length(tr_row),1)==repmat(tr_row',1,length(tr_row));        
        row_mat(eye(size(row_mat,1))==1)=0;
        col_mat(eye(size(col_mat,1))==1)=0;
        [~, row_ind] = find(max(row_mat));
        [~, col_ind] = find(max(col_mat));
        rm_vec = 1:length(tr_col);        
        dupe_indices = [dupe_indices tr_col(row_ind) tr_row(col_ind)];        
        tr_row = tr_row(~ismember(rm_vec,[row_ind col_ind]));
        tr_col = tr_col(~ismember(rm_vec,[row_ind col_ind]));                
    end    
    cat_fields = {'fluo','time','ap_vector','xPos','yPos'...
                    'cp_frames','all_frames'};
    for j = 1:length(tr_col)
        ID1 = min(tr_col(j),tr_row(j));
        ID2 = max(tr_col(j),tr_row(j));
        if ismember(ID1,remove_indices)            
            ID1 = match_indices(ID1==remove_indices);            
        end
        % take ID info from earlier trace (smaller ID)
        base = trace_struct(ID1);
        extra = trace_struct(ID2);        
        for f = 1:length(cat_fields)
            base.(cat_fields{f}) = [base.(cat_fields{f}) extra.(cat_fields{f})];
        end
        base.edge_flag = max(base.edge_flag,extra.edge_flag);        
        % Add particle info to nucleus struct
        fn = fieldnames(trace_struct);
        for f = 1:length(fn)
            trace_struct(ID1).(fn{f}) = base.(fn{f});
        end        
        remove_indices = [remove_indices ID2];
        match_indices = [match_indices ID1];
    end        
end

% remove extra entries
index_vector = 1:length(trace_struct);
nc_particles = [schnitz_struct.ParticleID];
tr_particles = [trace_struct.ParticleID];
trace_struct = trace_struct(~ismember(index_vector,[remove_indices dupe_indices]));
rm_particles = tr_particles([remove_indices dupe_indices]);
for i = 1:length(rm_particles)
    schnitz_struct(nc_particles==rm_particles(i)).ParticleID = NaN;
    schnitz_struct(i).xMean = mean(schnitz_struct(i).xPos);
    schnitz_struct(i).yMean = mean(schnitz_struct(i).yPos);
    schnitz_struct(i).xPosParticle = [];
    schnitz_struct(i).yPosParticle = [];
end
save([trace_name '.mat'],'trace_struct') 
save([nucleus_name '.mat'],'schnitz_struct') 

%% --------------Use 1D Clustering to Find Estimate Stripe Domains----------%
%%% color info
cm = jet(128);
increment = floor(size(cm,1)/7);
stripe_colors = cm(1+((1:7)-1)*increment,:);
cluster_sets = set_index;
    
%%% Use 1D Clustering to Find Stripe Centers
close all

stripe_radius = .015;
cluster_path = [ap_pos_path '/stripe_fits/'];
mkdir(cluster_path);


%%% tracking parameters
cluster_struct = struct;
t_window = 300;
track_times = 10:10:50;
% record param values
cluster_struct.track_times = track_times;
cluster_struct.t_window = t_window;
stripe_prior_mat = NaN(length(track_times),7);
%Create Arrays to store total fluorecence and # data points
%Careful with precision errors here...
ap_index = round(min([trace_struct.ap_vector]),2):.01:round(max([trace_struct.ap_vector]),2);
ap_index = round(ap_index,2);
ap_fluo_levels = zeros(length(track_times),length(ap_index));
%%% make smoothing kernel
sigma = 1;
sz = 2;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

%Iterate through sets
for j = 1:length(trace_struct)
    ap_path = trace_struct(j).ap_vector;
    fluo = trace_struct(j).fluo;
    time = trace_struct(j).time;
    time_filter = time(~isnan(fluo));
    fluo = fluo(~isnan(fluo));
    for t = 1:length(track_times)            
        start_time = track_times(t)*60-t_window;
        stop_time = track_times(t)*60+t_window;
        filter_vec = (time_filter >= start_time).*(time_filter < stop_time);        
        fluo_trunc = fluo(1==filter_vec); 
        ap_trunc = ap_path(1==filter_vec);
        for k = 1:length(fluo_trunc)                 
            ap = round(ap_trunc(k),2);
            %Add fluroescence value to appropriate bin
            ap_filter = ap_index==ap;            
            %Make sure that time point is assigned to 1 and only 1 AP
            %position
            if sum(ap_filter) ~= 1
                error('Error in AP Index Assignment');
            end
            ap_fluo_levels(t,ap_filter) = ap_fluo_levels(t,ap_filter) + fluo_trunc(k);                           
        end
    end
end

%%% Obtain Priors for stripe location using data distribution across sets
for t = 1:length(track_times)
    agg_fluo = ap_fluo_levels(t,:);
    %%% first tp
    agg_fluo_smooth = filter (gaussFilter,1, agg_fluo);
    raw_stripe_fig = figure;
    hold on
    plot(ap_index,ap_fluo_levels(t,:),'.')    
    y_lim = 1.1*max(ap_fluo_levels(t,:));
    plot(ap_index,agg_fluo_smooth,'LineWidth',2)
    grid on
    title(['Aggregated Fluorescence by AP (' num2str(track_times(t)) ')'])
    [s_ids, ~] = ginput;

    if track_times(t) < 20 % skip 5 for early periods
        s_ids = [s_ids(1:4); NaN ;s_ids(5:6)];
    end    
    for j = 1:length(s_ids)        
        plot([s_ids(j) s_ids(j)],[0 y_lim],'Color',stripe_colors(j,:),'LineWidth',1.5)
    end
    ylim([0 y_lim])    
    stripe_prior_mat(t,1:length(s_ids)) = s_ids';
    saveas(raw_stripe_fig, [ap_pos_path '/raw_stripes_t' num2str(track_times(t)) '.png'],'png')
    close all
end
cluster_struct.stripe_prior_mat = stripe_prior_mat;
cluster_struct.set_index = set_index;
%% perform clustering
% make smoothing kernel
sigma = 10;
sz = 10;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter);
%%% store results
centroid_1D_mat = NaN(length(track_times),7,length(set_index));
anterior_1D_mat = NaN(length(track_times),7,length(set_index));
posterior_1D_mat = NaN(length(track_times),7,length(set_index));
f_unit_cell = cell(length(track_times),length(set_index));
edge_flag_mat = NaN(length(track_times),7,length(set_index));
s_ids = 1:7;    
max_iters = 1000; %Max iterations allowed    
ap_size = .001; %Set granularity for ap stripe plots   
set_titles = {}; %Store set titles
set_subtitles = {};
%indicates flank
new_stripe_centers = NaN(length(set_index),7);        
% stripe_edge_flags = NaN(length(set_index),7);        
%Find stripe centers for each set
for s = 1:length(cluster_sets)
    set_fluo = [trace_struct(set_vec==cluster_sets(s)).fluo];
    set_ap = [trace_struct(set_vec==cluster_sets(s)).ap_vector];
    set_time = [trace_struct(set_vec==cluster_sets(s)).time];
    set_time = set_time(~isnan(set_fluo));
    set_fluo = set_fluo(~isnan(set_fluo));
    %set negatives to 0
    set_fluo(set_fluo < 0) = 0;
    %%% iterate through track times
    for t = 1:length(track_times)
        start_time = track_times(t)*60-t_window;
        stop_time = track_times(t)*60+t_window;
        t_filter = set_time>=start_time & set_time < stop_time;
        set_ap_t = set_ap(t_filter);
        set_fluo_t = set_fluo(t_filter);                 
        %Break cumulative fluo up into units of 10000. Each such unit is one
        %observation used for clustering (a way of coarse-grained, weighted
        %clustering)
        f_wt_vec = ceil(set_fluo_t/1e4);
        f_max = max(f_wt_vec);
        ap_wt_mat = repmat(set_ap_t,f_max,1);
        f_wt_mat = repmat(f_wt_vec,f_max,1);
        [~, ind_mat] = meshgrid(1:size(ap_wt_mat,2),1:size(ap_wt_mat,1));
        ap_wt_mat(f_wt_mat<ind_mat) = 0;
        f_unit_counts = reshape(ap_wt_mat,1,[]);
        f_unit_counts = f_unit_counts(f_unit_counts>0); 
        f_unit_cell{t,s} = f_unit_counts;
        %Cluster Data
        n_changes = 1;
        iter = 0;
        id_vec_old = zeros(1,length(f_unit_counts));                       
        %Remove stripe center priors that are not present in set
        new_stripe_centers = stripe_prior_mat(track_times==track_times(t),:);
        p_edge =  max(f_unit_counts);
        a_edge = min(f_unit_counts);
        stripe_id_vec = find((new_stripe_centers >= a_edge)...
                         &(new_stripe_centers <= p_edge));            
        new_stripe_centers = new_stripe_centers(stripe_id_vec);                        
        stripe_centers_orig = new_stripe_centers;            
        n_stripes = length(new_stripe_centers);
        %%% perform clustering
        while n_changes > 0                
            % find nearest center for each particle
            [m ,id_vec] = min(abs(repmat(f_unit_counts,length(new_stripe_centers),1)-new_stripe_centers'));                
            id_vec = stripe_id_vec(id_vec);            
            % check how many assignments changed
            n_changes = sum(id_vec~=id_vec_old);
            %Re-calculate centers
            stripe_centers_new = zeros(1,length(new_stripe_centers));
            for i = 1:length(new_stripe_centers)
                mean_p = nanmean(f_unit_counts(ismember(id_vec,stripe_id_vec(i))));
                stripe_centers_new(i) = mean_p;
            end                      
            iter = iter + 1;
            id_vec_old = id_vec;
            new_stripe_centers = stripe_centers_new;    
            if iter > max_iters
                warning('Maximum iterations exceeded');
                break
            end                
        end                            
        centroid_1D_mat(track_times==track_times(t),stripe_id_vec,s) = new_stripe_centers;          
%         stripe_edge_flags = new_stripe_centers<a_edge+stripe_radius|...
%                             new_stripe_centers>p_edge-stripe_radius;
%         edge_flag_mat(track_times==track_times(t),stripe_id_vec,s) = stripe_edge_flags;
%             centroid_1D_mat(track_times==track_times(t),stripe_edge_flags==1,s) = NaN; % remove edge cases
        %%% make fig
        fn = cp_filenames{s};
        tt_s_end = strfind(fn,'\Comp') - 1;
        tt_s_start = strfind(fn,'\2018') + 1;
        st_s_end = length(fn) - 4;
        t_string = fn(tt_s_start:tt_s_end);  
        st_string = fn(tt_s_end + 2: st_s_end);
        t_string = strrep(t_string,'_','-');
        st_string = strrep(st_string,'_','-');
        set_titles = [set_titles {t_string}];
        set_subtitles = [set_subtitles {st_string}];
        %Plot Stripes
        stripe_fig = figure('Visible', 'off');
        ap_vec = round(min(f_unit_counts),3):ap_size:round(max(f_unit_counts),3);
        hold on
        legend_string = {};
        b = [];
        for i = 1:n_stripes
            %plot full background                 
            hold on
            anterior_1D_mat(t,stripe_id_vec(i),s) = min(f_unit_counts(id_vec==stripe_id_vec(i)));
            posterior_1D_mat(t,stripe_id_vec(i),s) = max(f_unit_counts(id_vec==stripe_id_vec(i)));                
            ap_vec_plot_full = histc(f_unit_counts(id_vec==stripe_id_vec(i)),ap_vec);
            first_pt = find(ap_vec_plot_full,1);
            last_pt = find(ap_vec_plot_full,1,'last');
            ap_ct_plot = ap_vec(first_pt:last_pt);
            ap_vec_plot_full = ap_vec_plot_full(first_pt:last_pt);
            ap_vec_plot_full = conv(ap_vec_plot_full,gaussFilter,'same');
            norm_factor = conv(ones(size(ap_vec_plot_full)),gaussFilter,'same');
            ap_vec_plot_full = ap_vec_plot_full./norm_factor;

            b = [b bar(ap_ct_plot,ap_vec_plot_full,1,'FaceColor',...
                stripe_colors(stripe_id_vec(i),:),'FaceAlpha',.8,'EdgeAlpha',0)];     

            %overlay stripe centers
            legend_string = [legend_string {['Stripe ' num2str(stripe_id_vec(i))]}];               
        end
        grid on
        legend(b,legend_string{:});
        title(['Inferred Stripe Positions, ' t_string ' (Set ' num2str(s) ', Time ' num2str(track_times(t)) ')']);
        xlabel('AP Position (%)');
        ylabel('AU')                                 
        saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) 'time_' num2str(track_times(t)) '.png'],'png')
        saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) 'time_' num2str(track_times(t)) '.pdf'],'pdf')
        saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) 'time_' num2str(track_times(t)) '.fig'],'fig')
    end                     
end
%%% record boundary estimates
cluster_struct.anterior_1D_rough = anterior_1D_mat;
cluster_struct.posterior_1D_rough = posterior_1D_mat;
cluster_struct.centroid_1D_rough = centroid_1D_mat;
%%
%%% Review Stripe Classifications Manually
final_anterior_mat = anterior_1D_mat;
final_posterior_mat = posterior_1D_mat;
final_centroid_mat = centroid_1D_mat;
i = 1;
while i <= length(cp_filenames)    
    t = 1;
    while t <= length(track_times)
        close all
        f = openfig([cluster_path 'stripe_positions_set_' num2str(i) 'time_' num2str(track_times(t)) '.fig'],'visible');
        hold on
        [c_enter, ~] = ginput;        
        %if no corrections, proceed. Else implement changes and re-display    
        if ~isempty(c_enter)            
            new_stripe_id_vec = input('enter stripe IDs');
%             new_edge_flag_vec = input('flag edge cases');
            new_anterior_vec = NaN(1,7);
            new_posterior_vec = NaN(1,7);
            new_stripe_centers = NaN(1,7);
%             new_edge_flags = NaN(1,7);
            if length(new_stripe_id_vec)...
                    ~= length(c_enter)
                warning('inconsistent id and center vector sizes. retry')
                continue                
            end
            new_stripe_centers(ismember(1:7,new_stripe_id_vec)) = c_enter;
            
            f_unit_counts = f_unit_cell{t,i};                                    
            %%% make plot
            stripe_fig = figure;
            comp_centers = new_stripe_centers(~isnan(new_stripe_centers))';
            [m ,id_vec] = min(abs(repmat(f_unit_counts,length(comp_centers),1)-comp_centers));              
            stripe_vec_full = new_stripe_id_vec(id_vec);                
            %set x axis
            ap_vec = round(min(f_unit_counts),3):ap_size:round(max(f_unit_counts),3);
            hold on        
            legend_string = {};
            iter = 1;
            for k = new_stripe_id_vec
                legend_string = [legend_string {['Stripe ' num2str(k)]}];
                % record edges
                new_anterior_vec(k) = min(f_unit_counts(stripe_vec_full==k));
                new_posterior_vec(k) = max(f_unit_counts(stripe_vec_full==k));      
                ap_ct_plot = histc(f_unit_counts(stripe_vec_full==k),ap_vec);
                first_pt = find(ap_ct_plot,1);
                last_pt = find(ap_ct_plot,1,'last');
                ap_vec_plot = ap_vec(first_pt:last_pt);
                ap_ct_plot = ap_ct_plot(first_pt:last_pt);
                ap_ct_plot = conv(ap_ct_plot,gaussFilter,'same');
                norm_factor = conv(ones(size(ap_ct_plot)),gaussFilter,'same');
                ap_ct_plot = ap_ct_plot./norm_factor;

                bar(ap_vec_plot,ap_ct_plot,1,'FaceColor',...
                    stripe_colors(k,:),'FaceAlpha',.8,'EdgeAlpha',0);     
                               
                iter = iter + 1;
            end        
            grid on
            fn = cp_filenames{i};
            tt_s_end = strfind(fn,'\Comp') - 1;
            tt_s_start = strfind(fn,'\2017') + 1;
            st_s_end = length(fn) - 4;
            t_string = fn(tt_s_start:tt_s_end); 
            legend(legend_string{:});
            title(['Inferred Stripe Positions (Corrected), ' t_string ' (Set ' num2str(i) ' Time ' num2str(track_times(t)) ')']);
            xlabel('AP Position (%)');
            ylabel('AU')
            [x2, ~] = ginput;
            if isempty(x2)                
                final_centroid_mat(t,:,i) = new_stripe_centers;  
                final_anterior_mat(t,:,i) = new_anterior_vec;  
                final_posterior_mat(t,:,i) = new_posterior_vec;    
%                 final_edge_flag_mat(t,new_stripe_id_vec,i) = new_edge_flag_vec;
                saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(i) 'time_' num2str(track_times(t)) '_corrected.png'],'png')
                saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(i) 'time_' num2str(track_times(t)) '_corrected.pdf'],'pdf')
                saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(i) 'time_' num2str(track_times(t)) '_corrected.fig'],'fig')                          
                t = t + 1;
            end         
        else                    
            t = t + 1;
        end   
    end
    i = i + 1;
end
cluster_struct.final_centroid_mat = final_centroid_mat;
cluster_struct.final_anterior_mat = final_anterior_mat;
cluster_struct.final_posterior_mat = final_posterior_mat;
% cluster_struct.final_edge_flag_mat = final_edge_flag_mat;    
save([data_path 'stripe_clustering_results.mat'],'cluster_struct')

