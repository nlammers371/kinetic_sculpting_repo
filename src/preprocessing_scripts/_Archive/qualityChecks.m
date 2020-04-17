%------------------------Import Compiled Particles------------------------%
%Set path to folder containing relevant projects
folder_path = 'C:\Users\Nicholas\Dropbox (Garcia Lab)\mHMM\orig\';
% folder_path = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\DropboxSingleTraces\Eve2_ML';
project = 'mHMMeve2_orig_inf_set';
% outName = 'eve2Sets_2017_06_15_ml.mat'
outpath = [folder_path '/projects/' project '/' ];
% Keyword to ensure only sets from current project are pulled
keyword = '20sec';
%vector of data set numbers to include
include_vec = [9,10,19,20,21,22,23,24,26];
% exclude_vec = [12:20 5];
% exclude_vec(exclude_vec == 7) = 3;

%-------------------------Set Summary Parameters--------------------------%
nuclear_cycles = [14];
% Minimumt number of data points per summary stat
min_stat = 500;
%Define Relevant Grouping Regions
ap_grp_indices = {33:37,38:42,43:48};
ap_grp_names = {'Anterior Flank', 'Eve Stripe 2','Posterior Flank'};
 
if exist(outpath) ~= 7
    mkdir(outpath);
end
%NL: this could be made much better...
dir_struct = struct;
i_pass = 1;
dirinfo = dir(folder_path);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
subdirinfo = cell(length(dirinfo));
for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    % skip files lacking project keyword or containing names of skip sets
    if isempty(strfind(thisdir,keyword)) 
        continue
    end
    set_num_start_ind = strfind(thisdir,'_');
    set_num_start_ind = set_num_start_ind(end);
    set_num = str2num(thisdir(set_num_start_ind+1:end));
    
    if sum(set_num==include_vec) ~= 1 
        continue
    end
    
    particle_struct = dir(fullfile(folder_path,thisdir, 'CompiledParticles*'));
    ap_info_struct = dir(fullfile(folder_path,thisdir, 'APDetection'));
    dir_struct(i_pass).particles = particle_struct;
    dir_struct(i_pass).ap_info = ap_info_struct;
    dir_struct(i_pass).folder = thisdir;
    i_pass = i_pass + 1;
end

%struct to stor AP directories
ap_info_struct = struct;
filenames = {};
for i = 1:length(dir_struct)
    particle_struct = dir_struct(i).particles;
    if length(particle_struct) > 1
        warning(['Multiple Compile Particles Sets Detected for set ' num2str(include_vec(i))])
    end
    filenames = [filenames {[folder_path dir_struct(i).folder '\' particle_struct.name]}];
    ap_info = dir_struct(i).ap_info;
    for j = 1:length(ap_info)
        if strfind(ap_info(j).name,'HalfEmbryoArea') > 0
            ap_info_struct(i).SurfPath = [folder_path dir_struct(i).folder '\APDetection\' ap_info(j).name];
        elseif strfind(ap_info(j).name,'FullEmbryoArea') > 0           
            ap_info_struct(i).MidPath = [folder_path dir_struct(i).folder '\APDetection\' ap_info(j).name];
        end
    end
end
%%
%Data structure to store extracted trace sets
diff_traces = 0;
trace_struct = struct;
i_iter = 1;
for k = 1:length(filenames)
    raw_data  = load([filenames{k}]);
    time = raw_data.ElapsedTime*60;
%     Normalize to Start of 14
    
    traces_raw = raw_data.AllTracesVector;
    filter = sum(raw_data.ncFilter(:,ismember([raw_data.ncFilterID],nuclear_cycles)),2)>0;
    
    traces = traces_raw(:,filter==1);
    no_cycle_start = ~isnan(max(max(traces(time < 2*60,:)))) || max(max(traces(time < 2*60,:))) > 0;
    time = time - time(raw_data.nc14);
    for i = 1:size(traces,2)
        raw_trace = traces(:,i);
        if length(raw_trace(~isnan(raw_trace))) < 2
            continue
        end
        start = find(~isnan(raw_trace),1);
        stop = find(~isnan(raw_trace),1,'last');
        trunc_trace = [raw_trace(start:stop)'];
        trunc_time = time(start:stop); 
        [~, apPos] = max(raw_data.APFilter(i,:));
        trace_struct(i_iter).fluo = trunc_trace;
        trace_struct(i_iter).time = trunc_time;
        trace_struct(i_iter).AP = apPos;
        trace_struct(i_iter).set = filenames{k};
        trace_struct(i_iter).setID = k;
        trace_struct(i_iter).background = raw_data.MeanOffsetVector;
        trace_struct(i_iter).background_error = raw_data.SDOffsetVector;
        trace_struct(i_iter).background_time = time;
        trace_struct(i_iter).no_nc_start = no_cycle_start;
        i_iter = i_iter + 1;
    end
%     if k == 2
%         break
%     end
end
%% Histogram of Data Points by AP

sets = unique({trace_struct.set});
n_sets = length(sets);

set_titles = {};
for i = 1:length(sets)
    start = strfind(sets{i},'\2017') + 12;
    stop = strfind(sets{i},'\Compiled') - 1;
    string = sets{i};
    set_titles = {set_titles{:}  strrep(string(start:stop),'_',' ')};        
end

figure(1);
subplot(2,1,1);

colormap('jet');
cm = colormap;

% Unique AP values present in Data
ap_vec = unique([trace_struct.AP]);
%Store Trace and Data Point Counts By AP and Dataset
ap_ds_dp_mat = zeros(n_sets,length(ap_vec));
for a = 1:length(ap_vec)
    ap = ap_vec(a);
    for d = 1:n_sets
        ap_ds_dp_mat(d,a) = sum(length([trace_struct(1==([trace_struct.AP] == ap).*([trace_struct.setID] == d)).fluo]));        
    end
end
ap_ct_dp_vec = sum(ap_ds_dp_mat);

bar(ap_vec, ap_ct_dp_vec , 'Facecolor',cm(5,:),...
    'EdgeColor','black','FaceAlpha',.65,'BarWidth', 1);
title('Data Points by AP Region');
xlabel('AP Region (%)');
grid on
subplot(2,1,2);
histogram([trace_struct.AP],'Facecolor',cm(5,:));
title('Traces by AP Region');
xlabel('AP Region (%)');
grid on
saveas(figure(1), [outpath, 'data_histograms.png'],'png');

%%%Break out counts by Dataset
figure(2);
increment = floor(60 / n_sets);
% Array to store color mappings
set_colors = zeros(n_sets, 3);
hold on 
%Data Points
for i = 1:n_sets
    set_colors(i,:) = cm(1+(i-1)*increment,:);
    plot(ap_vec, ap_ds_dp_mat(i,:), '-o', 'Color',set_colors(i,:) ,'Linewidth', 1.5)
end
grid on
title('Data Points by AP Region and Set');
xlabel('AP Position (%)');
ylabel('# Points');
legend(set_titles{:});
hold off
saveas(figure(2), [outpath, 'data_count_plots.png'],'png');


%% Multi Hist Plots
close all
n_grps = length(ap_grp_indices);

increment = floor(60 / n_sets);
%Set size of fluo bins
granularity = 20;
max_fluo = ceil(max([trace_struct.fluo])/granularity)*granularity;
FluoBins = 1:granularity:max_fluo;
if exist([outpath 'fluo_his_plots/']) ~= 7
    mkdir([outpath 'fluo_his_plots/'])
end
fluo_his_fig = figure('Position',[0 0 1024 1024]);
% Struct to store hist infor for subsequent use
hist_info = struct;

for j = 1:n_sets
    for i = 1:n_grps
        ap_struct = trace_struct(ismember([trace_struct.AP],ap_grp_indices{i}));
        f_list = [];
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                fluo = ap_struct(a).fluo;
                f_list = [f_list fluo(fluo>0)];
            end
        end
        if isempty(f_list)
            continue
        end
        ap_ct = histc(f_list, FluoBins);        
        subplot(n_sets,n_grps, (j-1)*n_grps + i)
        b = bar(FluoBins, ap_ct / max(ap_ct), 'FaceColor',set_colors(j,:),...
            'EdgeColor',set_colors(j,:),'BarWidth', 1);
        set(gca,'fontsize',4)
        title(['Fluo, Set: ' set_titles{j} ' Region: ' ap_grp_names{i}]); %' Set:' sets{j}])
        axis([0,max_fluo,0 ,1])    
        grid on
        if i == 2;
            hist_info(j).hist_ct = ap_ct;
        end
    end
end
saveas(fluo_his_fig, [outpath, '_fluo_his.png'],'png');
hold off
%% Fig to Compare Embryo Orientation and Frame Position to Basic Trace Statistics
close all
if exist([outpath, 'Orientation'])~= 7
    mkdir([outpath, 'Orientation']);
end
%Get full distribution across all sets
eve2_fluo_list = [];
for i = 1:length(trace_struct);
    if ismember(trace_struct(i).AP, 38:42)
        eve2_fluo_list = [eve2_fluo_list trace_struct(i).fluo];
    end
end
        
for i = 1:n_sets
    image_fig = figure('Position', [0 0 1024 1024], 'Visible', 'off');
    subplot(2,2,1);
    imshow(ap_info_struct(i).MidPath)
    title(['Mid Saginal Alignment: Set ' set_titles{i}]);
    set(gca,'fontsize',6)
    
    subplot(2,2,3);
    imshow(ap_info_struct(i).SurfPath)
    title(['Surface Alignment: Set ' set_titles{i}]);
    set(gca,'fontsize',6)
    
    subplot(2,2,2);
    hold on
    full_dist = histc(eve2_fluo_list, FluoBins);
    bar(FluoBins, full_dist / sum(full_dist), 'FaceColor','black',...
            'EdgeColor','black', 'FaceAlpha', 0.5,'BarWidth', 1);
    set(gca,'fontsize',6)
    ap_ct = hist_info(i).hist_ct;
    bar(FluoBins, ap_ct / sum(ap_ct), 'FaceColor',set_colors(i,:),...
            'EdgeColor',set_colors(i,:),'FaceAlpha', 0.5, 'BarWidth', 1);
    set(gca,'fontsize',6)
%     title(['Fluo, Set: ' set_titles{j} ' Region: ' ap_grp_names{i}]); %' Set:' sets{j}])
    axis([0,max_fluo,0, max([full_dist / sum(full_dist) ap_ct / sum(ap_ct)])])  
    grid on
    title(['Fluorescent Intensities in Eve Stripe 2: Set ' set_titles{i}]);
    
    subplot(2,2,4);
    hold on
    all = area(ap_vec, sum(ap_ds_dp_mat,1)/sum(sum(ap_ds_dp_mat)));
    set(all,'FaceAlpha',0.5,'FaceColor','black');
    s = area(ap_vec, ap_ds_dp_mat(i,:)/sum(sum(ap_ds_dp_mat)));
    set(s,'FaceAlpha',0.5,'FaceColor',set_colors(i,:));
    grid on
    title(['Data Points by AP Position: Set ' set_titles{i}]);
    set(gca,'fontsize',10)
    saveas(image_fig, [outpath, 'Orientation/' 'orient_check_set_' num2str(include_vec(i)) '.png'],'png');
    hold off
end


% for i = 1:n_sets
%     
%     mid_fig = figure;
%     imshow(ap_info_struct(i).MidPath)
%     title(['Mid Saginal Alignment: Set ' set_titles{i}]);
%     set(gca,'fontsize',10)
%     saveas(mid_fig, [outpath, 'Orientation/' 'MidSag' num2str(include_vec(i)) '.png'],'png');
%     
%     surf_fig = figure;
%     imshow(ap_info_struct(i).SurfPath)
%     title(['Surface Alignment: Set ' set_titles{i}]);
%     set(gca,'fontsize',10)
%     saveas(surf_fig, [outpath, 'Orientation/' 'Surf_' num2str(include_vec(i)) '.png'],'png');
%     
%     f_fig = figure('Visible', 'off');
%     hold on
%     full_dist = histc(eve2_fluo_list, FluoBins);
%     bar(FluoBins, full_dist / sum(full_dist), 'FaceColor','black',...
%             'EdgeColor','black', 'FaceAlpha', 0.5,'BarWidth', 1);
%     set(gca,'fontsize',10)
%     ap_ct = hist_info(i).hist_ct;
%     bar(FluoBins, ap_ct / sum(ap_ct), 'FaceColor',set_colors(i,:),...
%             'EdgeColor',set_colors(i,:),'FaceAlpha', 0.5, 'BarWidth', 1);
%     set(gca,'fontsize',10)
% %     title(['Fluo, Set: ' set_titles{j} ' Region: ' ap_grp_names{i}]); %' Set:' sets{j}])
%     axis([0,max_fluo,0, max([full_dist / sum(full_dist) ap_ct / sum(ap_ct)])])  
%     grid on
%     title(['Fluorescent Intensities in Eve Stripe 2: Set ' set_titles{i}]);
%     saveas(f_fig, [outpath, 'Orientation/' 'Fluo Dist_' num2str(include_vec(i)) '.png'],'png');
%     hold off
%     
%     ap_fig = figure;
%     hold on
%     all = area(ap_vec, sum(ap_ds_dp_mat)/sum(sum(ap_ds_dp_mat)));
%     set(all,'FaceAlpha',0.5,'FaceColor','black');
%     s = area(ap_vec, ap_ds_dp_mat(i,:)/sum(sum(ap_ds_dp_mat)));
%     set(s,'FaceAlpha',0.5,'FaceColor',set_colors(i,:));
%     grid on
%     title(['Data Points by AP Position: Set ' set_titles{i}]);
%     set(gca,'fontsize',10)
%     saveas(ap_fig, [outpath, 'Orientation/' 'AP Points_' num2str(include_vec(i)) '.png'],'png');
% end

%% Multi HeatMap Plots
% Determine appropriate AP grouping
           
heat_fig = figure('Position',[0 0 2048 1024]);
for i = 1:n_grps
    ap_struct = trace_struct(ismember([trace_struct.AP],ap_grp_indices{i}));
    for j = 1:n_sets
        f_list = [];
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                fluo = ap_struct(a).fluo;
                f_list = [f_list fluo(fluo>0)];
            end
        end
        if isempty(f_list)
            continue
        end
        ap_ct = repmat(repelem(histc(f_list, FluoBins),granularity),1,1);        
        subplot(n_sets,n_grps, (j-1)*n_grps + i)
        sub_scheme = cm;
        colormap(sub_scheme);
       imagesc(ap_ct);
%         colorbar
        set(gca,'YTick',[])
        
        set(gca,'fontsize',6)
%         xlabel('Fluorescence (AU)')
        title(['Fluo, Set: '  set_titles{j} ' '  ap_grp_names{i}]); %' Set:' sets{j}])
%         axis([0,max_fluo,0 ,1])
    end
end
hold off
saveas(heat_fig, [outpath, 'fluo_heat.png'],'png');

%% Fluo Background Trends
close all
offset_fig = figure;
hold on
for i = 1:n_sets
    set_traces = trace_struct([trace_struct.setID]==i);
    plot(set_traces(1).background_time,set_traces(1).background, 'Color', set_colors(i,:),'Linewidth',1.5);
end
grid on
title ('Mean Offset Value over Time');
xlabel('Seconds (Relative to start of nc14)');
ylabel('AU (?)');
legend(set_titles{:});
saveas(offset_fig, [outpath, 'fluo_offset.png'],'png');
%% Cumulative Fluorescence within Eve Stripe 2 Region
close all
%Make Strings for legen entries
ap_range = ap_grp_indices{2};
figure(3)
ptile_list = zeros(1,n_sets);
% Set scale for x axis. Max of selected percentile across sets will be used
ptile = 97;
hold on
for j = 1:n_sets
    f_list = [];
    for a = 1:length(trace_struct)
        if strcmp(trace_struct(a).set,sets{j}) && ismember(trace_struct(a).AP,ap_range)
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
title(['Cumulative PDF (AP: ' num2str(min(ap_range)) '-' num2str(max(ap_range)) ')']); 
axis([0,max(ptile_list),0,1])
grid on
xlabel('AU')
legend(set_titles{:}, 'Location','southeast')
hold off
saveas(figure(3), [outpath, 'cum_his.png'],'png');

%% Plot Fluo Percentiles By AP and Data Set
med_fig = figure('Position',[0 0 1024 512]);

medians = nan(n_sets,length(ap_vec));
nTile40 = nan(n_sets,length(ap_vec));
nTile60 = nan(n_sets,length(ap_vec));
% hold on
for j = 1:n_sets
    for k = 1:length(ap_vec)
        f_list = [];
        ap_grp = ap_vec(k);
        ap_struct = trace_struct(ismember([trace_struct.AP],ap_grp));
        for a = 1:length(ap_struct)
            if strcmp(ap_struct(a).set,sets{j})
                f_list = [f_list ap_struct(a).fluo];
            end
        end
        if isempty(f_list)
            continue
        end
        f_list = f_list(f_list > 0);
%         f_list = f_list(f_list < 1500);
        medians(j,k) = median(f_list);
        nTile40(j,k) = prctile(f_list,40);
        nTile60(j,k) = prctile(f_list,60);
    end 
end

hold on
for i = 1:n_sets
    plot(ap_vec,medians(i,:), 'Color', set_colors(i,:),'Linewidth',2)
    legend(set_titles{:}, 'Location','southeast')
end
title('Median Fluorescent Intensities by AP Position and Data Set')
xlabel('AP Position (%)')
ylabel('AU');

for i = 1:n_sets
    p40_vec = nTile40(i,:); 
%     p40_vec = p40_vec(~isnan(p40_vec));
    p60_vec = nTile60(i,:);
%     p60_vec = p60_vec(~isnan(p60_vec));
%     x_vec = ap_vec(~isnan(p40_vec));         
    for j = 1:length(ap_vec)
        plot([ap_vec(j),ap_vec(j)], [p40_vec(j), p60_vec(j)], 'Color', cm(increment*(i-1)+1,:),'Linewidth',1)
    end
end

grid on
hold off
saveas(med_fig, [outpath, 'median_fluo_plots.png'],'png');

%% Mean Fluo By Region and Set
close all
for k = 1:length(ap_grp_indices)
    legend_names = {};
    fig =  figure('Position',[0 0 1024 512]);
    hold on
    ap_grp = ap_grp_indices{k};
    ap_struct = trace_struct(ismember([trace_struct.AP],ap_grp));
     
    for j = 1:n_sets
        set_traces = ap_struct([ap_struct.setID]==j);
        if isempty(set_traces)
            continue
        end
        if set_traces(1).no_nc_start == 1
            continue
        end
        legend_names = [legend_names set_titles{j}];
        set = set_traces(1).set;
        times = [set_traces.time];         
        fluo_values = [set_traces.fluo];
        fluo_values = fluo_values(~isnan(fluo_values));
        times = times(~isnan(fluo_values));
        unique_times = sort(unique(times));
        f_series = zeros(1,length(unique_times));
        for t = 1:length(unique_times)
            f_series(t) = mean(fluo_values(times == unique_times(t)));
        end
        plot(unique_times, f_series, 'Color', set_colors(j,:),'Linewidth',2);
    end
    title(['Mean Fluorescence by Data Set: ' ap_grp_names{k}]) %' (AP: ' num2str(ap_grp) ')'])
    xlabel('seconds');
    ylabel('AU');
    
    grid on;
    legend(legend_names{:}, 'Location','southwest')
    saveas(fig, [outpath, 'mean_temp_fluo_' ap_grp_names{k} '.png'],'png');
    
end



