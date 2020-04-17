% Script to generate pixel maps
% space and time
close all
% clear 

% set filenames
project = 'eve7stripes_inf_2018_04_28'; %Project Identifier
fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../../dat/' project '/']; % data mat directory

trace_name = [data_path 'raw_traces_01_' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_01_' project]; % names for compiled elipse struct
cluster_name = [data_path 'stripe_clustering_results.mat'];

fov_save_name = [data_path 'fov_partitions.mat'];
% load datasets
load(trace_name); % particle info
load(nucleus_name);
load(cluster_name);

xDim = 1024;
yDim = 256;
track_times = cluster_struct.track_times;
set_index = cluster_struct.set_index;
nc_set_vec = [schnitz_struct.setID];
% generate pixel-stripe maps
plot_times = min(track_times):max(track_times);
fov_stripe_partitions = struct;
ap_to_x_factors = NaN(1,length(set_index));
% set_vec = [trace_struct_final.setID];
for i = 1:length(set_index)    
    set_stripe_map = NaN(yDim,xDim,length(plot_times));
    set_centroids_pix = NaN(length(plot_times),7);
    set_centroids_ap = NaN(length(plot_times),7);
    
    % check for inversions
    all_ap = [trace_struct([trace_struct.setID]==i).ap_vector];
    all_x = [trace_struct([trace_struct.setID]==i).xPos];
    x_min = all_x(all_ap==min(all_ap));
    x_max = all_x(all_ap==max(all_ap));    
    
    flip_flag = 0;
    if x_min>x_max % deal with inversions
        flip_flag = 1;
%         all_x = xDim - all_x + 1;
    end   
    [x_ref, y_ref] = meshgrid(1:xDim,1:yDim);
    APAngle = trace_struct([trace_struct.setID]==i).APAngle;
    APLength = trace_struct([trace_struct.setID]==i).APLength;
    coordAZoom = trace_struct([trace_struct.setID]==i).coordAZoom;    
    Angle=atan2((y_ref-coordAZoom(2)),(x_ref-coordAZoom(1)));            
    Distance=sqrt((coordAZoom(2)-y_ref).^2+(coordAZoom(1)-x_ref).^2);
    APPosition=Distance.*cos(Angle-APAngle);
    APPosImage=APPosition/APLength;    
    ap_x_projection = mean(APPosImage);
    ap_x_pixel = length(ap_x_projection)/(max(ap_x_projection)-min(ap_x_projection)); % this only "works" becaus x axis is very 
                                                                                      % very nearly aligned with AP axis
    ap_x_pixel1 = APLength/100;
    a_mat_pix = NaN(size(cluster_struct.final_anterior_mat(:,:,i)));
    p_mat_pix = NaN(size(cluster_struct.final_anterior_mat(:,:,i)));
    c_mat_pix = NaN(size(cluster_struct.final_anterior_mat(:,:,i)));   
    for j = 1:size(a_mat_pix,1)
        for k = 1:size(a_mat_pix,2)
            [~, a_mat_pix(j,k)] = min(abs(ap_x_projection-cluster_struct.final_anterior_mat(j,k,i)));            
            [~, p_mat_pix(j,k)] = min(abs(ap_x_projection-cluster_struct.final_posterior_mat(j,k,i)));            
            [~, c_mat_pix(j,k)] = min(abs(ap_x_projection-cluster_struct.final_centroid_mat(j,k,i)));            
        end
    end        
%     error('afa')
    if flip_flag
        a_mat_pix = xDim - a_mat_pix + 1;
        p_mat_pix = xDim - p_mat_pix + 1;
        c_mat_pix = xDim - c_mat_pix + 1;
        APPosImage = fliplr(APPosImage);
    end
    c_mat_ap = cluster_struct.final_centroid_mat(:,:,i);    
    edge_flags = cluster_struct.final_edge_flag_mat(:,:,i);
    if i == 10
        edge_flags(:,4) = 1; % spot fix for set 10
    elseif i == 11
        edge_flags(:,4) = 0; % ditto for 11
    end
    % remove edge cases
    a_mat_pix(edge_flags==1|isnan(edge_flags)) = NaN;
    p_mat_pix(edge_flags==1|isnan(edge_flags)) = NaN;
    c_mat_pix(edge_flags==1|isnan(edge_flags)) = NaN;
    c_mat_ap(edge_flags==1|isnan(edge_flags)) = NaN;
    % store interp results
    a_mat_pix_interp = NaN(length(plot_times),7);
    p_mat_pix_interp = NaN(length(plot_times),7);
    c_mat_pix_interp = NaN(length(plot_times),7);        
    c_mat_ap_interp = NaN(length(plot_times),7);
    
    % interpolate
    for k = 1:7
        a_mat_pix_interp(:,k) = round(interp1(track_times,a_mat_pix(:,k)',plot_times));
        p_mat_pix_interp(:,k) = round(interp1(track_times,p_mat_pix(:,k),plot_times));
        c_mat_pix_interp(:,k) = round(interp1(track_times,c_mat_pix(:,k),plot_times));
        c_mat_ap_interp(:,k) = interp1(track_times,c_mat_ap(:,k),plot_times);        
    end    
    for t = 1:size(a_mat_pix_interp,1) % iterate through times
        a_vec = a_mat_pix_interp(t,:);
        p_vec = p_mat_pix_interp(t,:);
        c_vec_pix = c_mat_pix_interp(t,:);
        c_vec_ap = c_mat_ap_interp(t,:);
        stripe_id_vec = find(~isnan(c_vec_ap));
        a_vec = a_vec(~isnan(p_vec));        
        p_vec = p_vec(~isnan(p_vec));        
        c_vec_pix = c_vec_pix(~isnan(c_vec_pix));        
        
        set_centroids_pix(t,stripe_id_vec) = c_vec_pix;
        set_centroids_ap(t,:) = c_vec_ap;
        for j = 1:length(p_vec)            
            set_stripe_map(:,a_vec(j):p_vec(j),t) = stripe_id_vec(j);
        end           
    end       
    fov_stripe_partitions(i).ap_x_pixel = ap_x_pixel;
    fov_stripe_partitions(i).ap_x_pixel1 = ap_x_pixel1;
%     error('afs')
    fov_stripe_partitions(i).plot_times = plot_times;    
    fov_stripe_partitions(i).stripe_id_mat = set_stripe_map;
    fov_stripe_partitions(i).stripe_centroids = set_centroids_pix;
    fov_stripe_partitions(i).stripe_centroids_ap = set_centroids_ap;
    fov_stripe_partitions(i).ap_ref_mat = APPosImage;
    fov_stripe_partitions(i).t_track = track_times;
    fov_stripe_partitions(i).edge_flag_mat = edge_flags;
end
save([data_path 'fov_partitions.mat'],'fov_stripe_partitions')
save(nucleus_name, 'schnitz_struct');
%%
% %%%-------------------------QC Plots------------------------------------%%%
% set_titles = {};
% for i = 1:length(cp_filenames)
%     set_titles = [set_titles{:} {[' ' num2str(i) ' ']}];
% end
% %Make Color Palettes for use in figures
% precision = .01;
% n_sets = length(cp_filenames);
% cm = jet(128);
% increment = floor(size(cm,1) / n_sets);
% %Array to store color mappings
% set_colors = zeros(n_sets, 3);
% for i = 1:n_sets
%     set_colors(i,:) = cm(1+(i-1)*increment,:);
% end
% 
% %Set dimensions for figs 
% xDim = 1;
% yDim = 1;
% toggle = 1;
% while xDim*yDim < n_sets
%     if toggle == 1
%         xDim = xDim + 1;
%     else
%         yDim = yDim + 1;
%     end
%     toggle = toggle == 0;
% end
% 
% %Generate "coarse" AP vectors to aggregate fine-grained results
% %Have to be careful with floating point errors here
% coarse_ap_index = round(floor(round(ap_index,4)/precision+10e-6)*precision,2);
% coarse_ap = round(unique(coarse_ap_index),2);
% coarse_f_avg = zeros(size(ap_fluo_levels,1),length(coarse_ap));
% coarse_f_cum = zeros(size(ap_fluo_levels,1),length(coarse_ap));
% %Aggregate AP Fluo levels and TP Counts
% for i = 1:length(coarse_ap)
%    ap = coarse_ap(i);
%    coarse_filter = coarse_ap_index==ap;
%    if sum(coarse_filter) > precision / .0001 || sum(coarse_filter) < 1
%        error('Problem with AP Coarse Graining');
%    end
%    coarse_f_avg(:,i) = sum(ap_fluo_levels(:,coarse_ap_index==ap),2) ./ sum(ap_tp_counts(:,coarse_ap_index==ap),2);
%    coarse_f_cum(:,i) = sum(ap_fluo_levels(:,coarse_ap_index==ap),2);
% end
% 
% %-------------------------AP Averages-------------------------------------%
% mean_fluo_fig = figure('Position',[0 0 1536 1536]);
% % Struct to store hist infor for subsequent use
% max_mean = .9*max(coarse_f_avg(:));
% for j = 1:n_sets        
%     subplot(xDim,yDim,j)
%     hold on
%     bar(coarse_ap, nanmean(coarse_f_avg)  , 'FaceColor','black',...
%         'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);
%     bar(coarse_ap, coarse_f_avg(j,:)  , 'FaceColor',set_colors(j,:),...
%         'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);
%     set(gca,'fontsize',4)
%     title(strvcat(['Mean Fluorescence per Time Step NC 14, Set: ' set_titles{j}],...
%           set_subtitles{j})); %' Set:' sets{j}])
%     axis([min(coarse_ap),max(coarse_ap),0 ,max_mean])    
%     grid on
% end
% saveas(mean_fluo_fig, [ap_pos_path, 'mean_fluo_ap.png'],'png');
% 
% 
% % Integrated Fluorescence With Stripe Centers
% cumulative_fluo_fig = figure('Position',[0 0 1536 1536]);
% max_cum = max(coarse_f_cum(:));
% for j = 1:n_sets        
%     subplot(xDim,yDim, j)
%     hold on     
%     %Plot average profile
%     bar(coarse_ap, nanmean(coarse_f_cum)  , 'FaceColor','black',...
%         'FaceAlpha', .3,'EdgeColor','black','EdgeAlpha',.3,'BarWidth', 1);        
%     %Plot set profile
%     bar(coarse_ap, coarse_f_cum(j,:)  , 'FaceColor',set_colors(j,:),...
%         'FaceAlpha', .5,'EdgeColor',set_colors(j,:),'BarWidth', 1);        
%     set(gca,'fontsize',4)
%     title(strvcat(['Cumulative Fluorescence in NC 14, Set: ' set_titles{j}], ...
%         set_subtitles{j}),'Fontsize',6); %' Set:' sets{j}])
%     axis([min(coarse_ap),max(coarse_ap),0 ,max_cum])    
%     grid on
% end
% saveas(cumulative_fluo_fig, [ap_pos_path, 'cumulative_fluo_ap.png'],'png');
% 
% %--------------------------Multi Fluo Hist Plots--------------------------%
% % Set size of fluo bins
% fluopath = [fig_path '/fluo_statistics/'];
% if exist(fluopath) ~= 7
%     mkdir(fluopath)
% end
% max_fluo = ceil(max([trace_struct.fluo])); 
% granularity = floor(max_fluo/200);
% %Get set of unique stripes in trace_struct
% stripe_set = unique([trace_struct.stripe_id_coarse]);
% stripe_set = stripe_set(~isnan(stripe_set));
% for s = 1:length(stripe_set)
%     stripe = stripe_set(s);
%     
%     % Struct to store hist infor for subsequent use        
%     stripe_struct = trace_struct([trace_struct.stripe_id_coarse]==stripe);
%     s_sets = unique([stripe_struct.setID]);
%     stripe_titles = set_titles(s_sets);
%     stripe_subtitles = set_subtitles(s_sets);
%     stripe_colors = set_colors(s_sets,:);
%     n_sets_f = length(s_sets);
%     yDim_fluo = ceil(n_sets_f/xDim);
%     max_fluo = ceil(max([stripe_struct.fluo]));    
%     FluoBins = 1:granularity:ceil(max_fluo/granularity)*granularity;
%     fluo_his_fig = figure('Position',[0 0 xDim*256 yDim_fluo*256]);
%     for j = 1:n_sets_f      
%         f_list = [];
%         for a = 1:length(stripe_struct)
%             if stripe_struct(a).setID==s_sets(j)
%                 fluo = stripe_struct(a).fluo;
%                 f_list = [f_list fluo(fluo>0)];
%             end
%         end
%         if isempty(f_list)
%             continue
%         end
%         ap_ct = histc(f_list, FluoBins);        
%         subplot(yDim_fluo,xDim, j)
%         b = bar(FluoBins, ap_ct / max(ap_ct), 'FaceColor',stripe_colors(j,:),...
%             'EdgeColor',stripe_colors(j,:),'BarWidth', 1);
%         set(gca,'fontsize',4)
%         title(strvcat(['Fluo Distribution ,Stripe ' num2str(stripe) ' Set: ' stripe_titles{j}],...
%             stripe_subtitles{j})); 
%         axis([0,max_fluo,0 ,1])    
%         grid on
%     end
%     saveas(fluo_his_fig, [fluopath 'set_fluo_stripe_' num2str(stripe) '.png'],'png');
%     hold off
% end
% 
% %--------------------Cumulative Fluorescence------------------------------%
% %Make Strings for legen entries
% cf_fig = figure('Position',[0 0 1024 1024]);
% ptile_list = zeros(1,n_sets);
% % Set scale for x axis. Max of selected percentile across sets will be used
% ptile = 97;
% max_fluo = ceil(max([trace_struct.fluo]));    
% FluoBins = 1:granularity:ceil(max_fluo/granularity)*granularity;
% line_types = {'-','--','.-'};
% set_line_types = {};
% for i = 1:length(cp_filenames)
%     set_line_types = [set_line_types{:} {line_types{1+mod(i,3)}}];
% end    
% hold on
% for j = 1:n_sets
%     f_list = [];
%     for a = 1:length(trace_struct)
%         if trace_struct(a).setID == j %&& trace_struct(a).APbinID==0
%             f_list = [f_list trace_struct(a).fluo];
%         end
%     end
%     ptile_list(j) = prctile(f_list,ptile); 
%     if isempty(f_list)
%         continue
%     end
%     ap_ct = histc(f_list, FluoBins);
%     plot(FluoBins, cumsum(ap_ct) / sum(ap_ct), line_types{1+mod(j,3)},'Color',set_colors(j,:),'LineWidth', 2);
% end 
% title('Cumulative PDF');
% axis([0,max(ptile_list),0,1])
% grid on
% xlabel('AU')
% legend(set_titles{:}, 'Location','southeast')
% hold off
% saveas(cf_fig, [fluopath, '/set_cum_fluo.png'],'png');
% 
% 
% n_boots = 20;
% stripe_id_vec = [trace_struct.stripe_id_coarse];
% stripe_sub_id_vec = [trace_struct.stripe_sub_id_coarse];
% set_vec = [trace_struct.setID];
% % set_fluo_mean_array = NaN(length(include_vec),45);
% % set_fluo_se_array = NaN(length(include_vec),45);
% for i = 1:7
%     p = [];    
%     sets = unique(set_vec(stripe_id_vec==i&stripe_sub_id_vec==0));    
%     temporal_fig = figure('Position',[0 0 1024 512]);
%     hold on
%     time_vec = 1:50;
%     f_avg = zeros(1,length(time_vec));
%     legend_string = {};
%     for j = 1:length(sets)
%         SetID = sets(j);
%         set_struct_tr = trace_struct((stripe_id_vec==i)&(stripe_sub_id_vec==0)&(set_vec==SetID));
%         if length([set_struct_tr.fluo]) < 100 % skip poorly represented sets
%             continue
%         end
%         time_list = ceil([set_struct_tr.time]/60);
%         fluo_set = [set_struct_tr.fluo];
%         n_dp = length(fluo_set);
%         sample_vec = 1:n_dp;
%         f_array = NaN(n_boots,length(time_vec));
%         for n = 1:n_boots
%             s_ids = randsample(sample_vec,n_dp,true);
%             fluo_boot = fluo_set(s_ids);
%             time_boot = time_list(s_ids);
%             for t = time_vec
%                 f_array(n,t) = nanmean(fluo_boot(ismember(time_boot,t-2:t+2)));
%             end
%         end
%         sfm = nanmean(f_array);
%         sfe = nanstd(f_array);
% %         set_fluo_mean_array(j,:) = sfm;
% %         set_fluo_se_array(j,:) = sfe;
%         plot([time_vec ; time_vec],[sfm-sfe ; sfm + sfe],'Color',set_colors(sets(j),:),'LineWidth',1.5)
%         p = [p plot(time_vec, sfm, set_line_types{sets(j)},'Color',set_colors(sets(j),:),'LineWidth',1.5)];                
%         legend_string = {legend_string{:} ['Set ' num2str(sets(j))]};
%     end
%     title(['Average Fluorescence Over Time (5 min moving avg), Stripe ' num2str(i)]);
%     legend(p,legend_string{:},'Location','northwest')    
%     saveas(temporal_fig,[fluo_path 'time_trends_stripe_' num2str(i) '.png'],'png')
% end
% 
% %%% Make Temporal Heat Maps of Mean Activity For Each Set
% Tres = 3; % in minutes
% APres = 1; % percent AP
% t_vec = 0:3:48;
% for i = 1:length(cp_filenames)    
%     time_vec = [trace_struct([trace_struct.setID]==i).time];
%     fluo_vec = [trace_struct([trace_struct.setID]==i).fluo];
%     time_vec = time_vec(~isnan(fluo_vec));
%     fluo_vec = fluo_vec(~isnan(fluo_vec));
%     ap_vec = round([trace_struct([trace_struct.setID]==i).ap_vector]*100);
% %     set_ap_list = min(ap_vec)-1:max(ap_vec)+1;
%     ap_ref_vec = 10:90;
%     ap_mat = zeros(length(t_vec)-1,length(ap_ref_vec));
%     for j = 2:length(t_vec)
%         t1 = t_vec(j-1)*60;
%         t2 = t_vec(j)*60;
%         for a = 1:length(ap_ref_vec)
%             ap = ap_ref_vec(a);
%             ap_mat(j-1,a) = sum(fluo_vec(ap_vec==ap & time_vec < t2 & time_vec >= t1));
%         end
%     end
%     heat_fig = figure;
%     colormap(parula(128));
%     imagesc(ap_mat);
%     h = colorbar;
%     set(gca,'ytick',1:length(t_vec),'yticklabels',t_vec)
%     set(gca,'xtick',0:5:90,'xticklabels',10:5:90)
%     title(['Tracking Transcription Activity Over Time Set: ' num2str(i)])
%     xlabel('AP position')
%     ylabel('minutes into nc14')
%     saveas(heat_fig,[fluo_path 'temporal_heatmap_set' num2str(i) '.png'],'png')
% end
%         
% close all