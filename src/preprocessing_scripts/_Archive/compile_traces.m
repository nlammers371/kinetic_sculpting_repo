% Script to Compile Data Sets and Find Stripe Centers
close all
clear 
%------------------------Set Path Specs, ID Vars------------------------%
FolderPath = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\';
project = 'eve7stripes_inf_2018_03_27'; %Project Identifier
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
        trace_struct(j_pass).xPos = xPos;
        trace_struct(j_pass).yPos = yPos;
        trace_struct(j_pass).ap_vector = ap_positions;
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
set_index = [trace_struct.setID];
set_vec = unique(set_index);
% stitch together overlaps
remove_indices = [];
match_indices = [];
dupe_indices = [];
for s = 1:length(set_vec)
    set_struct = trace_struct(set_index==set_vec(s));     
    base_indices = find(set_index==set_vec(s));
    start_t = [];
    stop_t = [];
    start_x = [];
    stop_x = [];
    start_y = [];
    stop_y = [];
    for i = 1:length(set_struct)
        time = set_struct(i).time;
        xPos = set_struct(i).xPos;
        yPos = set_struct(i).yPos;
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

%%
%--------------Map total fluorecence to AP position for each set----------%
% exclude 1st 15 minutes since no mature stripes are present
start_time = 25*60;
stop_time = 60*60;
%Create Arrays to store total fluorecence and # data points
%Careful with rounding errors here...
ap_index = round(min([trace_struct.ap_vector]),3):.001:round(max([trace_struct.ap_vector]),3);
ap_index = round(ap_index,3);
ap_fluo_levels = zeros(length(cp_filenames),length(ap_index));
ap_tp_counts = zeros(length(cp_filenames),length(ap_index));

%Iterate through sets
for i = 1:length(cp_filenames)    
    traces = trace_struct([trace_struct.setID]==i); %Filter for relevant traces
    %For each time point in each trace, assign F(t) to AP(t)
    for j = 1:length(traces)
        ap_path = traces(j).ap_vector;
        fluo = traces(j).fluo;
        time = traces(j).time;
        time_filter = time(~isnan(fluo));
        fluo = fluo(~isnan(fluo));
        filter_vec = (time_filter >= start_time).*(time_filter < stop_time);        
        fluo_trunc = fluo(1==filter_vec); 
        ap_trunc = ap_path(1==filter_vec); 
        for k = 1:length(fluo_trunc)
            ap = round(ap_trunc(k),3);
            %Add fluroescence value to appropriate bin
            ap_filter = ap_index==ap;            
            %Make sure that time point is assigned to 1 and only 1 AP
            %position
            if sum(ap_filter) ~= 1
                error('Error in AP Index Assignment');
            end
            ap_fluo_levels(i,ap_filter) = ap_fluo_levels(i,ap_filter) + fluo_trunc(k);
            %Add tp count to appropriate AP bin
            ap_tp_counts(i,ap_filter) = ap_tp_counts(i,ap_filter) + 1;%/length(fluo_trunc);            
        end
    end
end

%%% Obtain Priors for stripe location using data distribution across sets

agg_counts = nansum(ap_tp_counts);
sigma = 5;
sz = 10;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
agg_counts_smooth = filter (gaussFilter,1, agg_counts);
raw_stripe_fig = figure;
hold on
plot(ap_index,sum(ap_tp_counts),'.')
hold on
plot(ap_index,agg_counts_smooth,'LineWidth',2)
grid on
title('Aggregated Fluorescence by AP')
[stripe_priors, ~] = ginput;
saveas(raw_stripe_fig, [ap_pos_path '/raw_stripes.png'],'png')
close all
%Calling these by eye at the moment
stripe_priors = stripe_priors';
if length(stripe_priors) ~= 7
    error('incorrect number of stripe priors')
end

%% Use 1D Clustering to Find Stripe Centers (should be generalized to 2D)-%
cm = jet(128);
stripe_radius = .015; %Impose stripe width of 3% AP    
increment = floor(size(cm,1)/7);
cluster_path = [ap_pos_path '/stripe_fits/'];
mkdir(cluster_path);
use_old_boundaries = input('Use previous stripe boundaries (1=yes, 0=no)?');
if use_old_boundaries == 1
    load([data_path 'stripe_boundaries_' project '.mat'])
    load([data_path 'stripe_id_cell_' project '.mat'])
else
    s_ids = 1:7;    
    max_iters = 1000; %Max iterations allowed    
    ap_size = .001; %Set granularity for ap stripe plots   
    set_titles = {}; %Store set titles
    set_subtitles = {};
    %indicates flank
    stripe_positions = cell(1,length(cp_filenames));
    stripe_id_cell = cell(1,length(cp_filenames)); 
    f_unit_cell = cell(1,length(cp_filenames)); 
    stripe_edge_cell = cell(1,length(cp_filenames)); 
    %Find stripe centers for each set
    for s = 1:length(cp_filenames)
        ap_array = ap_fluo_levels(s,:);
        %set negatives to 0
        ap_array(ap_array < 0) = 0;
        %find location of first and last nonzero elements
        start = find(ap_array,1);
        stop = find(ap_array,1,'last');
        bounded_ap = ap_index(start:stop);
        ap_array = ap_array(start:stop);
        f_unit_counts = [];
        %Break cumulative fluo up into units of 50000. Each such unit is one
        %observation used for clustering (a way of coarse-grained, weighted
        %clustering)
        for i = 1:length(ap_array)
            f = ceil(ap_array(i)/5e4);    
            f_unit_counts = [f_unit_counts repelem(bounded_ap(i),f)];   
        end
        
        %Cluster Data
        n_changes = 1;
        iter = 0;
        id_vec_old = zeros(1,length(f_unit_counts));
        %Remove stripe center priors that are not present in set
        stripe_centers = stripe_priors(1==((stripe_priors >= min(bounded_ap))...
            .*(stripe_priors <= max(bounded_ap))));
        %check for stripes on edge of FOV...let's not use these
        stripe_edge_flags = zeros(1,length(stripe_centers));
        stripe_edge_flags(stripe_centers-min(bounded_ap)<stripe_radius) = 1;
        stripe_edge_flags(-stripe_centers+max(bounded_ap)<stripe_radius) = 1;
        %identify stripes that appear to be in set
        stripe_id_vec = s_ids(ismember(stripe_priors,stripe_centers));
        n_stripes = length(stripe_centers);
        %perform clustering
        while n_changes > 0    
            id_vec = zeros(1,length(f_unit_counts));
            %for each particle, find nearest center
            for p = 1:length(f_unit_counts)
                ap = f_unit_counts(p);
                [m ,id] = nanmin(abs(ap-stripe_centers));
                id_vec(p) = id;
            end
            %check how many assignments changed
            n_changes = sum(id_vec~=id_vec_old);
            %Re-calculate centers
            stripe_centers_new = zeros(1,length(stripe_centers));
            for i = 1:length(stripe_centers)
                mean_p = nanmean(f_unit_counts(ismember(id_vec,i)));
                stripe_centers_new(i) = mean_p;
            end        
            iter = iter + 1;
            id_vec_old = id_vec;
            stripe_centers = stripe_centers_new;    
            if iter > max_iters
                warning('Maximum iterations exceeded');
                break
            end
        end
        %convenience vector
        index_vec = 1:length(id_vec);
        stripe_vec = zeros(1,length(id_vec));
        stripe_vec_full = zeros(1,length(id_vec));
        stripe_boundaries = [];
        %assign traces to stripe regions
        for i = 1:n_stripes
            sc = stripe_centers(i);
            distances = f_unit_counts(id_vec==i) - sc;
            
            if i == 1
                neg_edge = min(distances);
            else
                neg_edge = stripe_centers(i-1)+last_pos_edge - ...
                    stripe_centers(i);
            end
            pos_edge = max(distances);
            last_pos_edge = pos_edge;
            candidate_ids = index_vec(id_vec==i);        
            %ID AP positions that are in stripe center
            stripe_ids = candidate_ids(abs(distances) <= stripe_radius);
            stripe_vec(stripe_ids) = i;     
            %Allow radii of boundary regions to be flexible for time being
            stripe_vec_full(candidate_ids) = i;     
            stripe_boundaries = [stripe_boundaries; sc+neg_edge sc-stripe_radius sc sc+stripe_radius sc+pos_edge];       
        end
        
        
        stripe_positions{s} = stripe_boundaries(stripe_edge_flags==0,:);
        use_stripes = stripe_id_vec(stripe_edge_flags==0);
        stripe_id_cell{s} = stripe_id_vec(ismember(stripe_id_vec,use_stripes));
        f_unit_cell{s} = f_unit_counts;
        stripe_edge_cell{s} = stripe_edge_flags;
%         stripe_id_full_cell{s} = f_unit_counts;
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
            ap_vec_plot_full = histc(f_unit_counts(stripe_vec_full==i),ap_vec);        
            b = [b bar(ap_vec,ap_vec_plot_full,1,'FaceColor',cm(1+(stripe_id_vec(i)-1)*increment,:),'FaceAlpha',.2)];            
            %overlay stripe centers
            legend_string = [legend_string {['Stripe ' num2str(stripe_id_vec(i))]}];
            if stripe_edge_flags(i) == 1
                continue
            end
            ap_vec_plot = histc(f_unit_counts(stripe_vec==i),ap_vec);        
            bar(ap_vec,ap_vec_plot,1,'FaceColor',cm(1+(stripe_id_vec(i)-1)*increment,:),'FaceAlpha',.7)
        end
        grid on
        legend(b,legend_string{:});
        title(['Inferred Stripe Positions, ' t_string ' (Set ' num2str(s) ')']);
        xlabel('AP Position (%)');
        ylabel('AU')
        saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) '.png'],'png')
        saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) '.eps'],'epsc')
        saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) '.fig'],'fig')
    end
end


%% Review Stripe Classifications
i = 1;
while i <= length(cp_filenames)    
    close all
    f = openfig([cluster_path 'stripe_positions_set_' num2str(i) '.fig'],'visible');
    hold on
    [x, ~] = ginput;
    %if no corrections, proceed. Else implement changes and re-display    
    if ~isempty(x)
        f_unit_vec = f_unit_cell{i};
        new_stripe_id_vec = zeros(size(f_unit_vec));
        new_stripe_id_vec_full = zeros(size(f_unit_vec));
        new_stripe_ids = zeros(1,length(x));
        s_boundaries = zeros(length(x),5);
        s_boundaries(:,3) = x;
        s_boundaries(:,2) = x - stripe_radius;
        s_boundaries(:,4) = x + stripe_radius;
        s_boundaries(1,1) = x(1) - 2*stripe_radius;
        s_boundaries(end,end) = x(end) + 2*stripe_radius;
        for j = 1:length(x)
            new_center = x(j);
            [~, stripe_num] = min(abs(stripe_priors - new_center));
            new_stripe_ids(j) = stripe_num;
            if j ~= length(x)
                midpoint =x(j) + (x(j+1)-x(j))/2;
                s_boundaries(j,5) = midpoint;
                s_boundaries(j+1,1) = midpoint;
            end
            new_stripe_id_vec_full(1==((f_unit_vec < s_boundaries(j,end)).*...
                (f_unit_vec >= s_boundaries(j,1)))) = stripe_num;
            new_stripe_id_vec(1==((f_unit_vec < s_boundaries(j,4)).*...
                (f_unit_vec >= s_boundaries(j,2)))) = stripe_num;
        end
        stripe_fig = figure;
        %set x axis
        ap_vec = round(min(f_unit_vec),3):ap_size:round(max(f_unit_vec),3);
        hold on        
        legend_string = {};
        for k = new_stripe_ids
            ap_vec_plot_full = histc(f_unit_vec(new_stripe_id_vec_full==k),ap_vec);        
            p = bar(ap_vec,ap_vec_plot_full,1,'FaceColor',cm(1+(k-1)*increment,:),'FaceAlpha',.2);
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            legend_string = [legend_string {['Stripe ' num2str(k)]}];
            ap_vec_plot = histc(f_unit_vec(new_stripe_id_vec==k),ap_vec);        
            bar(ap_vec,ap_vec_plot,1,'FaceColor',cm(1+(k-1)*increment,:),'FaceAlpha',.7)
        end        
        grid on
        fn = cp_filenames{i};
        tt_s_end = strfind(fn,'\Comp') - 1;
        tt_s_start = strfind(fn,'\2017') + 1;
        st_s_end = length(fn) - 4;
        t_string = fn(tt_s_start:tt_s_end); 
        legend(legend_string{:});
        title(['Inferred Stripe Positions (Corrected), ' t_string ' (Set ' num2str(i) ')']);
        xlabel('AP Position (%)');
        ylabel('AU')
        [x2, ~] = ginput;
        if isempty(x2)
            stripe_positions{i} = s_boundaries;
            stripe_id_cell{i} = new_stripe_ids;
            saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(i) 'corrected.png'],'png')
            saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(i) '_corrected.eps'],'epsc')
            saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(i) '_corrected.fig'],'fig')           
            i = i + 1;
        end 
    else        
        i = i + 1;
    end   
end

%% ------------------------Assign Stripe Identity---------------------------%
for i = 1:length(trace_struct)
    setID = trace_struct(i).setID;
    s_boundaries = stripe_positions{setID};
    s_id_list = stripe_id_cell{setID};
    AP = mean(trace_struct(i).ap_vector);
    trace_struct(i).stripe_id_coarse = NaN;
    trace_struct(i).stripe_sub_id_coarse = NaN;
    trace_struct(i).stripe_radius = stripe_radius; 
    classified = 0;
    for j = 1:size(s_boundaries,1)
        sb = s_boundaries(j,:);
        if AP <= sb(5) && AP >= sb(1)
            trace_struct(i).stripe_id_coarse = s_id_list(j);            
            trace_struct(i).stripe_sub_id_coarse = 0;
            trace_struct(i).stripe_center_ap = sb(3);
            if AP >= sb(4)
                trace_struct(i).stripe_sub_id_coarse = 1;
            elseif AP <= sb(2)
                trace_struct(i).stripe_sub_id_coarse = -1;
            end
            classified = 1;
        end        
    end
end

% calculate ap to pixel calibrations for each set
pixel_to_ap_mat = NaN(length(cp_filenames),3);
for i = 1:length(cp_filenames)
    all_ap = [trace_struct([trace_struct.setID]==i).ap_vector];
    all_x = [trace_struct([trace_struct.setID]==i).xPos];
    all_y = [trace_struct([trace_struct.setID]==i).yPos];
    beta = regress(all_ap',[ones(size(all_x))' all_x' all_y']);
    pixel_to_ap_mat(i,:) = beta;
end

for i = 1:length(schnitz_struct)
    setID = schnitz_struct(i).setID;
    beta = pixel_to_ap_mat(setID,:);
    AP = beta(1) + beta(2)*schnitz_struct(i).xMean + beta(3)*schnitz_struct(i).yMean;    
    s_boundaries = stripe_positions{setID};
    s_id_list = stripe_id_cell{setID};    
    schnitz_struct(i).stripe_id_coarse = NaN;
    schnitz_struct(i).stripe_sub_id_coarse = NaN;
    schnitz_struct(i).stripe_radius = stripe_radius; 
    classified = 0;
    for j = 1:size(s_boundaries,1)
        sb = s_boundaries(j,:);
        if AP <= sb(5) && AP >= sb(1)
            schnitz_struct(i).stripe_id_coarse = s_id_list(j);            
            schnitz_struct(i).stripe_sub_id_coarse = 0;
            schnitz_struct(i).stripe_center_ap = sb(3);
            if AP >= sb(4)
                schnitz_struct(i).stripe_sub_id_coarse = 1;
            elseif AP <= sb(2)
                schnitz_struct(i).stripe_sub_id_coarse = -1;
            end
            classified = 1;
        end        
    end
end

%%% Make matrices assigning pixels to stripe region
xDim = 1024;
yDim = 256;
[x_ref, y_ref] = meshgrid(1:xDim,1:yDim);
fov_partitions = struct;
for i = 1:length(cp_filenames)
    beta = pixel_to_ap_mat(i,:);
    boundary_array = stripe_positions{i};
    stripe_id_vec = stripe_id_cell{i};    
    pixel_ap_id_mat = beta(1) + beta(2)*x_ref + beta(3)*y_ref;    
    pixel_stripe_id_mat = NaN(yDim,xDim);
    for j = 1:length(stripe_id_vec)
        filter = pixel_ap_id_mat >= boundary_array(j,1) & ...
                    pixel_ap_id_mat < boundary_array(j,end);
        pixel_stripe_id_mat(filter) = stripe_id_vec(j);
    end
    fov_partitions(i).pixel_stripe_id_mat = pixel_stripe_id_mat;
    fov_partitions(i).pixel_ap_id_mat = pixel_ap_id_mat;
end
%%% save data
save([trace_name '.mat'],'trace_struct') 
save([nucleus_name '.mat'],'schnitz_struct') 
save([data_path 'stripe_boundaries_' project '.mat'],'stripe_positions')
save([data_path 'stripe_id_cell_' project '.mat'],'stripe_id_cell')
save([data_path 'fov_partitions_' project '.mat'],'fov_partitions')

%%%-------------------------QC Plots------------------------------------%%%
set_titles = {};
for i = 1:length(cp_filenames)
    set_titles = [set_titles{:} {[' ' num2str(i) ' ']}];
end
%Make Color Palettes for use in figures
precision = .01;
n_sets = length(cp_filenames);
cm = jet(128);
increment = floor(size(cm,1) / n_sets);
%Array to store color mappings
set_colors = zeros(n_sets, 3);
for i = 1:n_sets
    set_colors(i,:) = cm(1+(i-1)*increment,:);
end

%Set dimensions for figs 
xDim = 1;
yDim = 1;
toggle = 1;
while xDim*yDim < n_sets
    if toggle == 1
        xDim = xDim + 1;
    else
        yDim = yDim + 1;
    end
    toggle = toggle == 0;
end

%Generate "coarse" AP vectors to aggregate fine-grained results
%Have to be careful with floating point errors here
coarse_ap_index = round(floor(round(ap_index,4)/precision+10e-6)*precision,2);
coarse_ap = round(unique(coarse_ap_index),2);
coarse_f_avg = zeros(size(ap_fluo_levels,1),length(coarse_ap));
coarse_f_cum = zeros(size(ap_fluo_levels,1),length(coarse_ap));
%Aggregate AP Fluo levels and TP Counts
for i = 1:length(coarse_ap)
   ap = coarse_ap(i);
   coarse_filter = coarse_ap_index==ap;
   if sum(coarse_filter) > precision / .0001 || sum(coarse_filter) < 1
       error('Problem with AP Coarse Graining');
   end
   coarse_f_avg(:,i) = sum(ap_fluo_levels(:,coarse_ap_index==ap),2) ./ sum(ap_tp_counts(:,coarse_ap_index==ap),2);
   coarse_f_cum(:,i) = sum(ap_fluo_levels(:,coarse_ap_index==ap),2);
end

%-------------------------AP Averages-------------------------------------%
mean_fluo_fig = figure('Position',[0 0 1536 1536]);
% Struct to store hist infor for subsequent use
max_mean = .9*max(coarse_f_avg(:));
for j = 1:n_sets        
    subplot(xDim,yDim,j)
    hold on
    bar(coarse_ap, nanmean(coarse_f_avg)  , 'FaceColor','black',...
        'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);
    bar(coarse_ap, coarse_f_avg(j,:)  , 'FaceColor',set_colors(j,:),...
        'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);
    set(gca,'fontsize',4)
    title(strvcat(['Mean Fluorescence per Time Step NC 14, Set: ' set_titles{j}],...
          set_subtitles{j})); %' Set:' sets{j}])
    axis([min(coarse_ap),max(coarse_ap),0 ,max_mean])    
    grid on
end
saveas(mean_fluo_fig, [ap_pos_path, 'mean_fluo_ap.png'],'png');


% Integrated Fluorescence With Stripe Centers
cumulative_fluo_fig = figure('Position',[0 0 1536 1536]);
max_cum = max(coarse_f_cum(:));
for j = 1:n_sets        
    subplot(xDim,yDim, j)
    hold on     
    %Plot average profile
    bar(coarse_ap, nanmean(coarse_f_cum)  , 'FaceColor','black',...
        'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);        
    %Plot set profile
    bar(coarse_ap, coarse_f_cum(j,:)  , 'FaceColor',set_colors(j,:),...
        'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);        
    set(gca,'fontsize',4)
    title(strvcat(['Cumulative Fluorescence in NC 14, Set: ' set_titles{j}], ...
        set_subtitles{j}),'Fontsize',6); %' Set:' sets{j}])
    axis([min(coarse_ap),max(coarse_ap),0 ,max_cum])    
    grid on
end
saveas(cumulative_fluo_fig, [ap_pos_path, 'cumulative_fluo_ap.png'],'png');

%--------------------------Multi Fluo Hist Plots--------------------------%
% Set size of fluo bins
fluopath = [fig_path '/fluo_statistics/'];
if exist(fluopath) ~= 7
    mkdir(fluopath)
end
max_fluo = ceil(max([trace_struct.fluo])); 
granularity = floor(max_fluo/200);
%Get set of unique stripes in trace_struct
stripe_set = unique([trace_struct.stripe_id_coarse]);
stripe_set = stripe_set(~isnan(stripe_set));
for s = 1:length(stripe_set)
    stripe = stripe_set(s);
    
    % Struct to store hist infor for subsequent use        
    stripe_struct = trace_struct([trace_struct.stripe_id_coarse]==stripe);
    s_sets = unique([stripe_struct.setID]);
    stripe_titles = set_titles(s_sets);
    stripe_subtitles = set_subtitles(s_sets);
    stripe_colors = set_colors(s_sets,:);
    n_sets_f = length(s_sets);
    yDim_fluo = ceil(n_sets_f/xDim);
    max_fluo = ceil(max([stripe_struct.fluo]));    
    FluoBins = 1:granularity:ceil(max_fluo/granularity)*granularity;
    fluo_his_fig = figure('Position',[0 0 xDim*256 yDim_fluo*256]);
    for j = 1:n_sets_f      
        f_list = [];
        for a = 1:length(stripe_struct)
            if stripe_struct(a).setID==s_sets(j)
                fluo = stripe_struct(a).fluo;
                f_list = [f_list fluo(fluo>0)];
            end
        end
        if isempty(f_list)
            continue
        end
        ap_ct = histc(f_list, FluoBins);        
        subplot(yDim_fluo,xDim, j)
        b = bar(FluoBins, ap_ct / max(ap_ct), 'FaceColor',stripe_colors(j,:),...
            'EdgeColor',stripe_colors(j,:),'BarWidth', 1);
        set(gca,'fontsize',4)
        title(strvcat(['Fluo Distribution ,Stripe ' num2str(stripe) ' Set: ' stripe_titles{j}],...
            stripe_subtitles{j})); 
        axis([0,max_fluo,0 ,1])    
        grid on
    end
    saveas(fluo_his_fig, [fluopath 'set_fluo_stripe_' num2str(stripe) '.png'],'png');
    hold off
end

%--------------------Cumulative Fluorescence------------------------------%
%Make Strings for legen entries
cf_fig = figure('Position',[0 0 1024 1024]);
ptile_list = zeros(1,n_sets);
% Set scale for x axis. Max of selected percentile across sets will be used
ptile = 97;
max_fluo = ceil(max([trace_struct.fluo]));    
FluoBins = 1:granularity:ceil(max_fluo/granularity)*granularity;
line_types = {'-','--','.-'};
set_line_types = {};
for i = 1:length(cp_filenames)
    set_line_types = [set_line_types{:} {line_types{1+mod(i,3)}}];
end    
hold on
for j = 1:n_sets
    f_list = [];
    for a = 1:length(trace_struct)
        if trace_struct(a).setID == j %&& trace_struct(a).APbinID==0
            f_list = [f_list trace_struct(a).fluo];
        end
    end
    ptile_list(j) = prctile(f_list,ptile); 
    if isempty(f_list)
        continue
    end
    ap_ct = histc(f_list, FluoBins);
    plot(FluoBins, cumsum(ap_ct) / sum(ap_ct), line_types{1+mod(j,3)},'Color',set_colors(j,:),'LineWidth', 2);
end 
title('Cumulative PDF');
axis([0,max(ptile_list),0,1])
grid on
xlabel('AU')
legend(set_titles{:}, 'Location','southeast')
hold off
saveas(cf_fig, [fluopath, '/set_cum_fluo.png'],'png');


n_boots = 20;
stripe_id_vec = [trace_struct.stripe_id_coarse];
stripe_sub_id_vec = [trace_struct.stripe_sub_id_coarse];
set_vec = [trace_struct.setID];
% set_fluo_mean_array = NaN(length(include_vec),45);
% set_fluo_se_array = NaN(length(include_vec),45);
for i = 1:7
    p = [];    
    sets = unique(set_vec(stripe_id_vec==i&stripe_sub_id_vec==0));    
    temporal_fig = figure('Position',[0 0 1024 512]);
    hold on
    time_vec = 1:50;
    f_avg = zeros(1,length(time_vec));
    legend_string = {};
    for j = 1:length(sets)
        SetID = sets(j);
        set_struct = trace_struct((stripe_id_vec==i)&(stripe_sub_id_vec==0)&(set_vec==SetID));
        if length([set_struct.fluo]) < 100 % skip poorly represented sets
            continue
        end
        time_list = ceil([set_struct.time]/60);
        fluo_set = [set_struct.fluo];
        n_dp = length(fluo_set);
        sample_vec = 1:n_dp;
        f_array = NaN(n_boots,length(time_vec));
        for n = 1:n_boots
            s_ids = randsample(sample_vec,n_dp,true);
            fluo_boot = fluo_set(s_ids);
            time_boot = time_list(s_ids);
            for t = time_vec
                f_array(n,t) = nanmean(fluo_boot(ismember(time_boot,t-2:t+2)));
            end
        end
        sfm = nanmean(f_array);
        sfe = nanstd(f_array);
%         set_fluo_mean_array(j,:) = sfm;
%         set_fluo_se_array(j,:) = sfe;
        plot([time_vec ; time_vec],[sfm-sfe ; sfm + sfe],'Color',set_colors(sets(j),:),'LineWidth',1.5)
        p = [p plot(time_vec, sfm, set_line_types{sets(j)},'Color',set_colors(sets(j),:),'LineWidth',1.5)];                
        legend_string = {legend_string{:} ['Set ' num2str(sets(j))]};
    end
    title(['Average Fluorescence Over Time (5 min moving avg), Stripe ' num2str(i)]);
    legend(p,legend_string{:},'Location','northwest')    
    saveas(temporal_fig,[fluo_path 'time_trends_stripe_' num2str(i) '.png'],'png')
end

%%% Make Temporal Heat Maps of Mean Activity For Each Set
Tres = 3; % in minutes
APres = 1; % percent AP
t_vec = 0:3:48;
for i = 1:length(cp_filenames)    
    time_vec = [trace_struct([trace_struct.setID]==i).time];
    fluo_vec = [trace_struct([trace_struct.setID]==i).fluo];
    time_vec = time_vec(~isnan(fluo_vec));
    fluo_vec = fluo_vec(~isnan(fluo_vec));
    ap_vec = round([trace_struct([trace_struct.setID]==i).ap_vector]*100);
%     set_ap_list = min(ap_vec)-1:max(ap_vec)+1;
    ap_ref_vec = 10:90;
    ap_mat = zeros(length(t_vec)-1,length(ap_ref_vec));
    for j = 2:length(t_vec)
        t1 = t_vec(j-1)*60;
        t2 = t_vec(j)*60;
        for a = 1:length(ap_ref_vec)
            ap = ap_ref_vec(a);
            ap_mat(j-1,a) = sum(fluo_vec(ap_vec==ap & time_vec < t2 & time_vec >= t1));
        end
    end
    heat_fig = figure;
    colormap(parula(128));
    imagesc(ap_mat);
    h = colorbar;
    set(gca,'ytick',1:length(t_vec),'yticklabels',t_vec)
    set(gca,'xtick',0:5:90,'xticklabels',10:5:90)
    title(['Tracking Transcription Activity Over Time Set: ' num2str(i)])
    xlabel('AP position')
    ylabel('minutes into nc14')
    saveas(heat_fig,[fluo_path 'temporal_heatmap_set' num2str(i) '.png'],'png')
end
        
