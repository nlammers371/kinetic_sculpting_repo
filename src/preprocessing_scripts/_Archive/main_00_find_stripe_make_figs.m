%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
% folder_path = 'C:\Users\Nicholas\Dropbox (Garcia Lab)\eve2spots\';
% folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\eve2spots\';
folder_path = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\';
%Assign Project Identifier
project = 'eve7stripes_inf_2018_02_20';

figpath = ['../../fig/experimental_system/preprocessing/' project '/'];
datapath = ['../../dat/' project '/'];
if exist(datapath) ~= 7
    mkdir(datapath);
end
if exist(figpath) ~= 7
    mkdir(figpath);
end
ap_pos_path = [figpath 'ap_positioning/'];
if exist(ap_pos_path) ~= 7
    mkdir(ap_pos_path)
end  
% Keyword to ensure only sets from current project are pulled
keyword = '_30uW_550V';
%vector of data set numbers to include

%-------------------------Set Summary Parameters--------------------------%
%Pull CompiledParticle Sets from all relevant project folders
dir_struct = struct;
i_pass = 1;
dirinfo = dir(folder_path);
%remove non-directories
dirinfo(~[dirinfo.isdir]) = [];  
subdirinfo = cell(length(dirinfo));
for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    % skip files lacking project keyword or containing names of skip sets
    if isempty(strfind(thisdir,keyword)) 
        continue
    end    
    particle_struct = dir(fullfile(folder_path,thisdir, 'CompiledParticles*'));    
    dir_struct(i_pass).particles = particle_struct;
    dir_struct(i_pass).folder = thisdir;
    i_pass = i_pass + 1;
end
%%
%Save Filepaths of all CompiledParticle Sets
filenames = {};
APnames = {};
for i = 1:length(dir_struct)
    particle_struct = dir_struct(i).particles;
    if isempty(particle_struct)
        continue    
    elseif length(particle_struct) > 1
        warning(['Multiple Compile Particles Sets Detected for set '...
            num2str(i) ', taking first'])
    end    
    filenames = [filenames {[folder_path dir_struct(i).folder ...
            '\' particle_struct(1).name]}];
    APnames = [APnames {[folder_path dir_struct(i).folder '\APDetection.mat']}];    
end

n_sets = length(filenames);
%Data structure to store extracted trace sets
trace_struct = struct;
j_pass = 1;
for i = 1:n_sets
    %AP Info
    load(APnames{i})
    %Angle between the x-axis and the AP-axis
    APAngle = round(atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)))*360 / (2*pi));    
    %Correction for if APAngle is in quadrants II or III
    if coordPZoom(1)-coordAZoom(1) < 0
        APAngle = APAngle + pi;
    end
    raw_data = load([filenames{i}]);
    time_raw = raw_data.ElapsedTime*60;
    particles = raw_data.CompiledParticles;
    traces_raw = raw_data.AllTracesVector;
    frames_raw = 1:length(time_raw);
    %Get frame that marks start of nc14
    start_frame = raw_data.nc14;
    %Filter trace, frame,and time arrays
    traces = traces_raw(start_frame:end,:);
    time = time_raw(start_frame:end);
    time = time - min(time);
    frames = frames_raw(start_frame:end);
    fn = filenames{i};
    %Iterate through traces
    for j = 1:size(traces,2)
        raw_trace = traces(:,j);     
        trace_start = find(~isnan(raw_trace),1);
        trace_stop = find(~isnan(raw_trace),1,'last');
        trace_full = raw_trace(trace_start:trace_stop)';
        trunc_trace = raw_trace(~isnan(raw_trace));
        %Remove dots
        if length(trunc_trace) < 2
            continue
        end        
        %remove NaN values
        trunc_time = time(~isnan(raw_trace));
        time_full = time(trace_start:trace_stop);
        trunc_frames = frames(~isnan(raw_trace));
        %Why the hell is there an APPos field with negative values?
        ap_positions_raw = particles(j).APpos;
        fov_xPos_raw = particles(j).xPos;
        fov_yPos_raw = particles(j).yPos;
        pt_frames = particles(j).Frame;
        ap_positions = ap_positions_raw(ismember(pt_frames,trunc_frames));        
        xPos = fov_xPos_raw(ismember(pt_frames,trunc_frames));
        yPos = fov_yPos_raw(ismember(pt_frames,trunc_frames));        
        trace_struct(j_pass).APAngle = APAngle;
        trace_struct(j_pass).xPos = xPos;        
        trace_struct(j_pass).yPos = yPos;        
        trace_struct(j_pass).ap_vector = ap_positions;  
        trace_struct(j_pass).MeanAP = mean(ap_positions);  
        trace_struct(j_pass).fluo = trunc_trace';
        trace_struct(j_pass).time = trunc_time;
        trace_struct(j_pass).fluo_full = trace_full;
        trace_struct(j_pass).time_full = time_full;
        trace_struct(j_pass).setID = i;
        trace_struct(j_pass).set = filenames{i};
        trace_struct(j_pass).FluoError = particles(j).FluoError;
        trace_struct(j_pass).OriginalParticle = particles(j).OriginalParticle;
        j_pass = j_pass + 1;
    end
end
%%
%--------------Map total fluorecence to AP position for each set----------%
%By Default Let's Assume that we wish to take full nc14
start_time = 0*60;
stop_time = 1000*60;
%Create Arrays to store total fluorecence and # data points
%Careful with rounding errors here...
ap_index = round(min([trace_struct.ap_vector]),4):.0001:round(max([trace_struct.ap_vector]),4);
ap_index = round(ap_index,4);
ap_fluo_levels = zeros(length(filenames),length(ap_index));
ap_tp_counts = zeros(length(filenames),length(ap_index));

%Iterate through sets
for i = 1:length(filenames)
    %Filter for relevant traces
    traces = trace_struct([trace_struct.setID]==i);
    %For each time point in each trace, assign F(t) to AP(t)
    for j = 1:length(traces)
        ap_path = traces(j).ap_vector;
        fluo = traces(j).fluo;
        time = traces(j).time;
        filter_vec = (time >= start_time).*(time < stop_time);
        fluo_trunc = fluo(1==filter_vec); 
        ap_trunc = ap_path(1==filter_vec); 
        for k = 1:length(fluo_trunc)
            ap = round(ap_trunc(k),4);
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
%% Obtain Priors for stripe location using data distribution across sets

colormap('jet');
cm = colormap;
increment = floor(size(cm,1)/7);

agg_fluo_levels = sum(ap_fluo_levels);
sigma = 50;
sz = 100;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize
agg_fluo_smooth = filter (gaussFilter,1, agg_fluo_levels);
raw_stripe_fig = figure;
hold on
plot(ap_index,sum(ap_fluo_levels),'.')
hold on
plot(ap_index,agg_fluo_smooth,'LineWidth',2)
grid on
title('Aggregated Fluorescence by AP')
[stripe_priors, ~] = ginput;
saveas(raw_stripe_fig, [ap_pos_path '/raw_stripes.png'],'png')
close all
%Calling these by eye at the moment
stripe_priors = stripe_priors';
% stripe_priors = [.32,.39,.47,.55,.63,.69,.75];
%% Use 1D Clustering to Find Stripe Centers (should be generalized to 2D)-%
cluster_path = [ap_pos_path '/stripe_fits/'];
if exist(cluster_path) ~= 7
    mkdir(cluster_path);
end
use_old_boundaries = input('Use previous stripe boundaries (1=yes, 0=no)?');
if use_old_boundaries == 1
    load([datapath 'stripe_boundaries_' project '.mat'])
    load([datapath 'stripe_id_cell_' project '.mat'])
else
    s_ids = 1:7;
    %Impose stripe width of 3 AP
    stripe_radius = .015;
    %Max iterations allowed
    max_iters = 1000;    
    %Set granularity for ap stripe plots
    ap_size = .001;
    %Store set titles
    set_titles = {};
    set_subtitles = {};
    %indicates flank)
    stripe_positions = cell(1,length(filenames));
    stripe_id_cell = cell(1,length(filenames)); 
    f_unit_cell = cell(1,length(filenames)); 
    %Find stripe centers for each set
    for s = 1:length(filenames)
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
            neg_edge = min(distances);
            pos_edge = max(distances);
            candidate_ids = index_vec(id_vec==i);        
            %ID AP positions that are in stripe center
            stripe_ids = candidate_ids(abs(distances) <= stripe_radius);
            stripe_vec(stripe_ids) = i;     
            %Allow radii of boundary regions to be flexible for time being
            stripe_vec_full(candidate_ids) = i;     
            stripe_boundaries = [stripe_boundaries; sc+neg_edge sc-stripe_radius sc sc+stripe_radius sc+pos_edge];       
        end
        stripe_positions{s} = stripe_boundaries;
        stripe_id_cell{s} = stripe_id_vec;
        f_unit_cell{s} = f_unit_counts;
    %     stripe_id_full_cell{s} = f_unit_counts;
        fn = filenames{s};
        tt_s_end = strfind(fn,'\Comp') - 1;
        tt_s_start = strfind(fn,'\2017') + 1;
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
        for i = 1:n_stripes
            %plot full background
            ap_vec_plot_full = histc(f_unit_counts(stripe_vec_full==i),ap_vec);        
            p = bar(ap_vec,ap_vec_plot_full,'FaceColor',cm(1+(stripe_id_vec(i)-1)*increment,:),'FaceAlpha',.2);
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            %overlay stripe centers
            legend_string = [legend_string {['Stripe ' num2str(stripe_id_vec(i))]}];
            ap_vec_plot = histc(f_unit_counts(stripe_vec==i),ap_vec);        
            bar(ap_vec,ap_vec_plot,'FaceColor',cm(1+(stripe_id_vec(i)-1)*increment,:),'FaceAlpha',.7)
        end
        grid on
        legend(legend_string{:});
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
while i <= length(filenames)    
    f = openfig([cluster_path 'stripe_positions_set_' num2str(i) '.fig'],'visible');
    hold on
    [x ~] = ginput;
    %if no corrections, proceed. Else implement changes and re-display
    close all
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
            p = bar(ap_vec,ap_vec_plot_full,'FaceColor',cm(1+(k-1)*increment,:),'FaceAlpha',.2);
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            legend_string = [legend_string {['Stripe ' num2str(k)]}];
            ap_vec_plot = histc(f_unit_vec(new_stripe_id_vec==k),ap_vec);        
            bar(ap_vec,ap_vec_plot,'FaceColor',cm(1+(k-1)*increment,:),'FaceAlpha',.7)
        end        
        grid on
        fn = filenames{i};
        tt_s_end = strfind(fn,'\Comp') - 1;
        tt_s_start = strfind(fn,'\2017') + 1;
        st_s_end = length(fn) - 4;
        t_string = fn(tt_s_start:tt_s_end); 
        legend(legend_string{:});
        title(['Inferred Stripe Positions (Corrected), ' t_string ' (Set ' num2str(i) ')']);
        xlabel('AP Position (%)');
        ylabel('AU')
        [x2 ~] = ginput;
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
%%
%------------------------Assign Stripe Identity---------------------------%
for i = 1:length(trace_struct)
    setID = trace_struct(i).setID;
    s_boundaries = stripe_positions{setID};
    s_id_list = stripe_id_cell{setID};
    AP = trace_struct(i).MeanAP;
    trace_struct(i).stripe_id = NaN;
    trace_struct(i).stripe_radius = stripe_radius;   
    for j = 1:size(s_boundaries,1)
        sb = s_boundaries(j,:);
        if AP <= sb(5) && AP >= sb(1)
            trace_struct(i).stripe_id = s_id_list(j);            
            trace_struct(i).stripe_sub_id = 0;
            trace_struct(i).stripe_center_ap = sb(3);
            if AP >= sb(4)
                trace_struct(i).stripe_sub_id = 1;
            elseif AP <= sb(2)
                trace_struct(i).stripe_sub_id = -1;
            end
        end        
    end           
end
%Remove nans (temporary)
trace_struct = trace_struct(~isnan([trace_struct.stripe_id]));
if length(trace_struct(isnan([trace_struct.stripe_id]))) > 5
    error('>5 nans in stripe assignment');
end

save([datapath 'raw_traces_' project '.mat'],'trace_struct') 
save([datapath 'stripe_boundaries_' project '.mat'],'stripe_boundaries')
save([datapath 'stripe_id_cell_' project '.mat'],'stripe_id_cell')
%-------------Plot Mean and Cumulative Fluorescence by Set----------------%
%Make Titles for Plots (This will lonly work for eve2 format set titles
%currently)
 
%Set dimensions for figs (Need to automate this)
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
% Set granularity (.01 = 1% AP precision)
precision = 0.01;

%Make Color Palettes for use in figures
colormap('jet')
cm = colormap;
increment = floor(size(cm,1) / n_sets);
%Array to store color mappings
set_colors = zeros(n_sets, 3);
for i = 1:n_sets
    set_colors(i,:) = cm(1+(i-1)*increment,:);
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
max_mean = max(max(coarse_f_avg));
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
max_cum = 1.1*max(max(coarse_f_cum));
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
fluopath = [figpath '/fluo_statistics/'];
if exist(fluopath) ~= 7
    mkdir(fluopath)
end
max_fluo = ceil(max([trace_struct.fluo])); 
granularity = floor(max_fluo/200);
%Get set of unique stripes in trace_struct
stripe_set = unique([trace_struct.stripe_id]);
stripe_set = stripe_set(~isnan(stripe_set));
for s = 1:length(stripe_set)
    stripe = stripe_set(s);
    
    % Struct to store hist infor for subsequent use        
    stripe_struct = trace_struct([trace_struct.stripe_id]==stripe);
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
