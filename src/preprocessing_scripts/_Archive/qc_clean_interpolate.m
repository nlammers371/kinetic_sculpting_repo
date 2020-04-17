
%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
% folder_path = 'C:\Users\Nicholas\Dropbox (Garcia Lab)\eve2spots\';
% folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\eve2spots\';
folder_path = 'D:\Data\Augusto\LivemRNA\Data\Dropbox\eveProject\eve7stripes\';
%Assign Project Identifier
project = 'eve7stripes_inf';

figpath = ['../../fig/' project '/'];
datapath = ['../../dat/' project '/'];
ap_pos_path = [figpath 'ap_positioning/'];
if exist(ap_pos_path) ~= 7
    mkdir(ap_pos_path)
end  
% Keyword to ensure only sets from current project are pulled
keyword = 'eve';
%vector of data set numbers to include

%-------------------------Set Summary Parameters--------------------------%
if exist(datapath) ~= 7
    mkdir(datapath);
end
if exist(figpath) ~= 7
    mkdir(figpath);
end
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

%Save Filepaths of all CompiledParticle Sets
filenames = {};
APnames = {};
for i = 1:length(dir_struct)
    particle_struct = dir_struct(i).particles;
    if length(particle_struct) > 1
        warning(['Multiple Compile Particles Sets Detected for set ' num2str(i) '. Taking first instance'])
    end
    filenames = [filenames {[folder_path dir_struct(i).folder '\' particle_struct(1).name]}];
    APnames = [APnames {[folder_path dir_struct(i).folder '\APDetection.mat']}];
end
%%
n_sets = length(filenames);
%Data structure to store extracted trace sets
trace_struct = struct;
j_pass = 1;
for i = 1:length(filenames)
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
        j_pass = j_pass + 1;
    end
end
%% Compile Aggregate Fluroescence and Obs Counter Per AP
unique_ap = unique(round([trace_struct.ap_vector],2));
ap_sums = zeros(1,length(unique_ap));
ap_counts = zeros(1,length(unique_ap));
for i = 1:length(trace_struct)
    trace = trace_struct(i).fluo;
    ap = round(trace_struct(i).ap_vector,2);
    for j = 1:length(trace)
        ap_sums(unique_ap==ap(j)) = ap_sums(unique_ap==ap(j)) + trace(j);
        ap_counts(unique_ap==ap(j)) = ap_counts(unique_ap==ap(j)) + 1;
    end
end


%% Map total fluorecence to AP position for each set
%By Default Let's Assume that we wish to take full nc14
start_time = 0*60;
stop_time = 60*60;
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
%% Use 1D Clustering to Find Stripe Centers (should be generalized to 2D)
cluster_path = [ap_pos_path '/stripe_fits/'];
if exist(cluster_path) ~= 7
    mkdir(cluster_path);
end
%NL: I made these priors up. Should check/refine
stripe_priors = [.32,.39,.48,.56,.62,.68,.74];
s_ids = 1:7;
%impose maximum stripe radius (width = 2*r)
max_radius = .025;
%Max iterations allowed
max_iters = 1000;    
%Percentile cutoff for stripe centers and flanks
threshold = .68;
%Set granularity for ap stripe plots
ap_size = .001;
%Store set titles
set_titles = {};
%indicates flank)
stripe_positions = cell(1,length(filenames));
stripe_id_cell = cell(1,length(filenames)); 
%Find stripe centers for each set
for s = 1:length(filenames)
    ap_array = ap_fluo_levels(s,:);
    start = find(ap_array,1);
    stop = find(ap_array,1,'last');
    bounded_ap = ap_index(start:stop);
    ap_array = ap_array(start:stop);
    f_unit_counts = [];
    for i = 1:length(ap_array)
        f = ceil(ap_array(i)/10000);    
        f_unit_counts = [f_unit_counts repelem(bounded_ap(i),f)];   
    end
    %Cluster Data
    n_changes = 1;
    iter = 0;
    id_vec_old = zeros(1,length(f_unit_counts));
    stripe_centers = bounded_ap(ismember(bounded_ap,stripe_priors));
    stripe_id_vec = s_ids(ismember(stripe_centers,stripe_priors));
    n_stripes = length(stripe_centers);

    while n_changes > 0    
        id_vec = zeros(1,length(f_unit_counts));
        for p = 1:length(f_unit_counts)
            ap = f_unit_counts(p);
            [m ,id] = nanmin(abs(ap-stripe_centers));
            id_vec(p) = id;
        end
        n_changes = sum(id_vec~=id_vec_old);
        %Re-calculate centers
        stripe_centers_new = zeros(1,length(stripe_centers));
        for i = 1:length(stripe_centers)
            mean_p = nanmean(f_unit_counts(ismember(id_vec,i)));
            stripe_centers_new(i) = mean_p;
        end
        if iter > max_iters
            warning('Maximum iterations exceeded');
            break
        end
        iter = iter + 1;
        id_vec_old = id_vec;
        stripe_centers = stripe_centers_new;    
    end
    %convenience vector
    index_vec = 1:length(id_vec);
    stripe_vec = zeros(1,length(id_vec));
    stripe_boundaries = [];
    for i = 1:n_stripes
        sc = stripe_centers(i);
        distances = abs(f_unit_counts(id_vec==i) - sc);
        candidate_ids = index_vec(id_vec==i);
        sr = min(prctile(distances,threshold*100),max_radius);
        stripe_ids = candidate_ids(distances <= sr);
        stripe_vec(stripe_ids) = i;     
        stripe_boundaries = [stripe_boundaries; sc-sr sc sc + sr];       
    end
    stripe_positions{s} = stripe_boundaries;
    stripe_id_cell{s} = stripe_id_vec;
    fn = filenames{s};
    tt_s_end = strfind(fn,'\Comp') - 1;
    tt_s_start = strfind(fn,'\2017') + 1;
    t_string = fn(tt_s_start:tt_s_end);  
    set_titles = [set_titles {t_string}];
    %Plot Stripes
    stripe_fig = figure('Visible', 'off');
    colormap('jet');
    cm = colormap;
    increment = floor(size(cm,1)/(n_stripes+1));
    ap_vec = round(min(f_unit_counts),3):ap_size:round(max(f_unit_counts),3);
    hold on
    for i = 0:n_stripes
        ap_vec_plot = histc(f_unit_counts(stripe_vec==i),ap_vec);        
        bar(ap_vec,ap_vec_plot,'FaceColor',cm(1+(1+i-1)*increment,:),'FaceAlpha',.5)
    end
    grid on
    title(['Inferred Stripe Positions, ' t_string ' (Set ' num2str(s) ')']);
    xlabel('AP Position (%)');
    ylabel('AU')
    saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) '.png'],'png')
    saveas(stripe_fig, [cluster_path 'stripe_positions_set_' num2str(s) '.eps'],'epsc')
end
%% Assign Stripe Identity

for i = 1:length(trace_struct)
    setID = trace_struct(i).setID;
    s_boundaries = stripe_positions{setID};
    s_id_list = stripe_id_cell{setID};
    AP = trace_struct(i).MeanAP;
    trace_struct(i).stripe_id = NaN;
    trace_struct(i).stripe_radius = NaN;
    trace_struct(i).stripe_threshold = threshold;
    for j = 1:size(s_boundaries,1)
        sb = s_boundaries(j,:);
        if AP <= sb(3) && AP >= sb(1)
            trace_struct(i).stripe_id = s_id_list(j);
            trace_struct(i).stripe_radius = sb(2) - sb(1);
        end        
    end       
end
save([datapath 'raw_traces_' project '.mat'],'trace_struct') 

%% Plot Mean and Cumulative Fluorescence by Set
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
mean_fluo_fig = figure('Position',[0 0 1024 1024]);
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
    title(['Mean Fluorescence per Time Step NC 14, Set: ' set_titles{j}]); %' Set:' sets{j}])
    axis([min(coarse_ap),max(coarse_ap),0 ,max_mean])    
    grid on
end
saveas(mean_fluo_fig, [ap_pos_path, 'mean_fluo_ap.png'],'png');

% Integrated Fluorescence With Stripe Centers
cumulative_fluo_fig = figure('Position',[0 0 1024 1024]);
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
    title(['Cumulative Fluorescence in NC 14, Set: ' set_titles{j}],'Fontsize',6); %' Set:' sets{j}])
    axis([min(coarse_ap),max(coarse_ap),0 ,max_cum])    
    grid on
end
saveas(cumulative_fluo_fig, [ap_pos_path, 'cumulative_fluo_ap.png'],'png');

% %%
% %Mean Fluo Over Time
% %Structure to store total fluo and fluo counts for each AP and Set
% mean_fluo_struct = struct;
% max_time = floor(max([trace_struct.time]/60));
% min_ap = 34;
% max_ap = 45;
% bin_size = 2;
% all_ap = (min_ap+.5):2:(max_ap-.5);
% aggregate_ap_counts = zeros(max_time+1,length(all_ap));
% aggregate_ap_fluo = zeros(max_time+1,length(all_ap));
% for k = 1:length(filenames)
%     set_traces = trace_struct([trace_struct.setID]==k);    
%     ap_counts = zeros(max_time+1,length(all_ap));
%     ap_fluo = zeros(max_time+1,length(all_ap));
%     for j = 1:length(set_traces)       
%         t = floor(set_traces(j).time/60);
%         f = set_traces(j).fluo;
%         ap = floor(100*set_traces(j).ap_vector);
%         for i = 1:length(f)
%             if ismember(ap(i),floor(all_ap)) || ismember(ap(i),ceil(all_ap)) 
%                 ap_counts(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap)))...
%                     = ap_counts(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap))) + 1;
%                 ap_fluo(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap))) ...
%                     = ap_fluo(t(i)+1,1==(ap(i)==floor(all_ap))+(ap(i)==ceil(all_ap))) + f(i);
%             end
%         end
%     end
%     aggregate_ap_counts = aggregate_ap_counts + ap_counts;
%     aggregate_ap_fluo = aggregate_ap_fluo + ap_fluo;
%     mean_fluo_struct(k).ap_vec = all_ap;
%     mean_fluo_struct(k).ap_counts = ap_counts;
%     mean_fluo_struct(k).ap_fluo = ap_fluo;
%     mean_fluo_struct(k).set_id = k;
% end
% % Get average overal behavior
% aggregate_fluo_mean = sum(aggregate_ap_fluo,2) ./ sum(aggregate_ap_counts,2);
% 
% % Make Fig
% close all
% temporal_fluo_fig = figure('Position',[0 0 2048 2048]);
% n_sets = length(filenames);
% 
% increment = floor(60/length(all_ap));
% index_vector = 1 + increment*(all_ap - min(all_ap))/2;
% colormap('winter');
% cm = colormap;
% 
% legend_string = {};
% 
% for i = 1:n_sets
%     subplot(dims,dims,i);
%     hold on
%     mean_fluo = mean_fluo_struct(i).ap_fluo ./ mean_fluo_struct(i).ap_counts;
%     for j = 1:size(mean_fluo,2)
%         plot(mean_fluo(:,j),'Color',[cm(index_vector(all_ap(j)==all_ap),:) .5],'Linewidth',1.5);            
%     end
%     plot(aggregate_fluo_mean, '-', 'Color', 'black','Linewidth',1.5)
%     axis([0 50 0 1000])
%     grid on
%     xlabel('Minutes into nc14');
%     ylabel('Mean Fluorescence (AU)');
%     title(['Mean Fluorescence Over Time by AP: Set ' num2str(i)],'FontSize',8);
%     set(gca,'FontSize',6)
% end
% 
% for ap = all_ap
%     ap_string = ['AP ' num2str(ap)];    
%     legend_string = {legend_string{:} ap_string};
% end
% legend(legend_string{:});
% 
% saveas(temporal_fluo_fig, [ap_pos_path, 'mean_fluo_temp.png'],'png');
%% Multi Fluo Hist Plots
close all
%Set size of fluo bins
max_fluo = ceil(max([trace_struct.fluo]));
granularity = floor(max_fluo/200);
FluoBins = 1:granularity:ceil(max_fluo/granularity)*granularity;
fluopath = [figpath '/fluo_statistics/'];
if exist(fluopath) ~= 7
    mkdir(fluopath)
end
fluo_his_fig = figure('Position',[0 0 1024 1024]);
% Struct to store hist infor for subsequent use
hist_info = struct;
% Just plot all AP bins for each set for now
for j = 1:n_sets
    bin_struct = trace_struct;%trace_struct([trace_struct.APbinID]==0);
    f_list = [];
    for a = 1:length(bin_struct)
        if bin_struct(a).setID==j
            fluo = bin_struct(a).fluo;
            f_list = [f_list fluo(fluo>0)];
        end
    end
    if isempty(f_list)
        continue
    end
    ap_ct = histc(f_list, FluoBins);        
    subplot(xDim,yDim, j)
    b = bar(FluoBins, ap_ct / max(ap_ct), 'FaceColor',set_colors(j,:),...
        'EdgeColor',set_colors(j,:),'BarWidth', 1);
    set(gca,'fontsize',4)
    title(['Fluo Distribution in Stripe Center, Set: ' set_titles{j}]); 
    axis([0,max_fluo,0 ,1])    
    grid on
end

saveas(fluo_his_fig, [fluopath 'set_stripe_fluo.png'],'png');
hold off

%% Cumulative Fluorescence within Eve Stripe 2 Region
close all
%Make Strings for legen entries
cf_fig = figure;
ptile_list = zeros(1,n_sets);
% Set scale for x axis. Max of selected percentile across sets will be used
ptile = 97;
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
    plot(FluoBins, cumsum(ap_ct) / sum(ap_ct), 'Color',set_colors(j,:),'LineWidth', 2);
end 
title('Cumulative PDF for Stripe Centers');
axis([0,max(ptile_list),0,1])
grid on
xlabel('AU')
legend(set_titles{:}, 'Location','southeast')
hold off
saveas(cf_fig, [fluopath, '/set_cum_fluo.png'],'png');
