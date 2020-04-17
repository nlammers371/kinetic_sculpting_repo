% Script to track stripe activity centers and (appx) mRNA levels across
% space and time
close all
clear 

% set filenames
project = 'eve7stripes_inf_2018_03_27'; %Project Identifier

fig_path = ['../../fig/experimental_system/' project '/preprocessing/'];
data_path = ['../../dat/' project '/']; % data mat directory

flux_dynamics_path = [fig_path 'ap_positioning/stripe_dynamics/'];
mkdir(flux_dynamics_path);
trace_name = [data_path 'raw_traces_' project]; % names for compiled trace struct
nucleus_name = [data_path 'ellipse_info_' project]; % names for compiled elipse struct
fov_name = [data_path 'fov_partitions_' project '.mat'];
stripe_save_name = [data_path 'stripe_pos_' project '.mat'];
% load datasets
load(trace_name); % particle info
load(fov_name); % ap and stripe info at pixel level
load(nucleus_name);
save_spline_figs = 1; % if 1 generates spline figs
%% generate fluorescence maps 
%%% Color Info
cm = jet(128);
increment = floor(size(cm,1)/7);
%%% fit variables
max_disp = 60; % max permissible fluo dispersion for stripe classfication 
stripe_radius = 35; % pixels .015;
kernel_radius = 60; % radius of gauss kernel...nucleus diameter ~= 20-25
kernel_sigma = 25; % this is kind of arbitrary
[x_ref, y_ref] = meshgrid(1:2*kernel_radius+1,1:2*kernel_radius+1);
x_ref = x_ref - kernel_radius - 1;
y_ref = y_ref - kernel_radius - 1;
r_mat = sqrt(x_ref.^2 + y_ref.^2);
g_kernel = exp(-(r_mat/(2*kernel_sigma))); % gauss kernel

%%% time stuff
t_window = 2; % lag/lead averaging window size in minutes
t_vec = 25:50; %start fitting stripes at 25 min

%%% make indexing vectors for traces
xPos_vec_particle = [];
yPos_vec_particle = [];
mean_xPos_vec_particle = [];
mean_yPos_vec_particle = [];
set_vec_particle = [];
fluo_vec_particle = [];
time_vec_particle = [];
particle_vec = [];
for i = 1:length(trace_struct)
    fluo = trace_struct(i).fluo;
    time = trace_struct(i).time;
    time = time(~isnan(fluo));
    fluo = fluo(~isnan(fluo));    
    if isempty(fluo)
        error('asfa')
    end
    particle_vec = [particle_vec repelem(trace_struct(i).ParticleID,length(trace_struct(i).yPos))];
    fluo_vec_particle = [fluo_vec_particle fluo];
    time_vec_particle = [time_vec_particle time];
    xPos_vec_particle = [xPos_vec_particle trace_struct(i).xPos];
    yPos_vec_particle = [yPos_vec_particle trace_struct(i).yPos];
    mean_xPos_vec_particle = [mean_xPos_vec_particle mean(trace_struct(i).xPos(time>=600))];
    mean_yPos_vec_particle = [mean_yPos_vec_particle mean(trace_struct(i).yPos(time>=600))];
    set_vec_particle = [set_vec_particle repelem(trace_struct(i).setID, length(trace_struct(i).yPos))];    
end
% same for nuclei
xPos_vec_nc = [];
yPos_vec_nc = [];
mean_xPos_vec_nc = [];
mean_yPos_vec_nc = [];
set_vec_nc = [];
fluo_vec_nc = [];
time_vec_nc = [];
nc_vec = [];
for i = 1:length(schnitz_struct)    
    time = schnitz_struct(i).time;                
    nc_vec = [nc_vec repelem(schnitz_struct(i).ncID,length(schnitz_struct(i).yPos))];
    time_vec_nc = [time_vec_nc time];
    xPos_vec_nc = [xPos_vec_nc schnitz_struct(i).xPos];
    yPos_vec_nc = [yPos_vec_nc schnitz_struct(i).yPos];
    mean_xPos_vec_nc= [mean_xPos_vec_nc mean(schnitz_struct(i).xPos(time>=600))];
    mean_yPos_vec_nc = [mean_yPos_vec_nc mean(schnitz_struct(i).yPos(time>=600))];
    set_vec_nc = [set_vec_nc repelem(schnitz_struct(i).setID, length(schnitz_struct(i).yPos))];    
end

set_index = unique(set_vec_particle);
stripe_class_vec = NaN(1,length(trace_struct));
stripe_pos_struct = struct; % save stripe location arrays
for i = 1:length(set_index)
    xp_set_vec = xPos_vec_particle(set_vec_particle==set_index(i));
    yp_set_vec = yPos_vec_particle(set_vec_particle==set_index(i));
    fluo_set_vec = fluo_vec_particle(set_vec_particle==set_index(i));
    time_set_vec = time_vec_particle(set_vec_particle==set_index(i));
    ap_mat = fov_partitions(i).pixel_ap_id_mat;
    stripe_mat = fov_partitions(i).pixel_stripe_id_mat;
    %%% handle inversions
    mean_ap_vec = nanmean(ap_mat);
    if mean_ap_vec(1) > mean_ap_vec(end)
        stripe_mat = fliplr(stripe_mat);
        ap_mat = fliplr(ap_mat);
        xp_set_vec = size(stripe_mat,2) - xp_set_vec + 1;
        mean_xPos_vec_particle = size(stripe_mat,2) - mean_xPos_vec_particle + 1;
    end
    stripe_id_vec = unique(reshape(stripe_mat(~isnan(stripe_mat)),1,[]),'stable');
    temp_fluo_array = NaN(size(stripe_mat,1),size(stripe_mat,2),length(t_vec));
    temp_spline_mat = NaN(size(stripe_mat,1),length(stripe_id_vec),length(t_vec));
    temp_disp_mat = NaN(size(stripe_mat,1),length(stripe_id_vec),length(t_vec));
    for j = 1:length(t_vec)
        frame_array = zeros(size(stripe_mat));
        t = t_vec(j)*60;
        t_start = t - t_window*60;
        t_stop = t + t_window*60;
        t_fluo = fluo_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_x = xp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        t_y = yp_set_vec(time_set_vec>=t_start&time_set_vec<t_stop);
        idx = sub2ind(size(frame_array), t_y, t_x);
        frame_array(idx) = t_fluo;
        % apply smoothing kernel
        norm_ref_array = ones(size(frame_array));
        norm_ref_array = conv2(norm_ref_array,g_kernel,'same');
        gauss_array = conv2(frame_array,g_kernel,'same');
        gauss_array = gauss_array./norm_ref_array;
        centroid_mat = NaN(size(stripe_mat,1),length(stripe_id_vec)); % true centroids
        dispersion_mat = NaN(size(stripe_mat,1),length(stripe_id_vec)); % true centroids
        spline_mat = NaN(size(stripe_mat,1),length(stripe_id_vec)); % spline fits
        for k = 1:length(stripe_id_vec)
            stripe_id = stripe_id_vec(k);            
            x_vec = [];
            y_vec = [];
            for y = 1:size(stripe_mat,1)
                y_strip =  gauss_array(y,:);
                y_strip(stripe_mat(y,:)~=stripe_id) = 0;                
                indices = 1:size(stripe_mat,2);
                mi = sum(indices.*y_strip)/sum(y_strip);                
                centroid_mat(y,k) = mi;
                dispersion_mat(y,k) = std(indices,y_strip); % weighted std
                %%% apply weights
                wt_vec = floor(y_strip/100);
                if max(wt_vec) == 0
                    continue
                end
                ind_wt_mat = repmat(indices,max(wt_vec),1);
                wt_mat = repmat(wt_vec,max(wt_vec),1);
                [x_grid,y_grid] = meshgrid(1:size(wt_mat,2),1:size(wt_mat,1));
                ind_wt_mat(y_grid>wt_mat) = 0;
                
                ind_wt_mat = ind_wt_mat(ind_wt_mat>0);
                ind_wt_vec = reshape(ind_wt_mat,1,[]);
                x_vec = [x_vec ind_wt_vec];
                y_vec = [y_vec repelem(y,length(ind_wt_vec))];                
            end            
            center_vec = centroid_mat(:,k);
            y_vec_sp = (1:size(stripe_mat,1))';
%             f = fit(y_vec_sp(~isnan(center_vec)),centroid_mat(~isnan(center_vec),k),'smoothingspline',...
%                 'SmoothingParam',0.01);    
            pp = polyfit(y_vec,x_vec,4);
            poly = polyval(pp,unique(y_vec));
%             error('afsa')
            spline_mat(:,k) = poly;
        end
        if save_spline_figs
            centroid_fig = figure;
            centroid_fig.Visible = 'off';
            centroid_fig.Position = [100 100 1024 256];
            hold on
            imagesc(gauss_array);
            scatter(reshape(centroid_mat,1,[]),repmat(1:size(stripe_mat,1),1,length(stripe_id_vec)),...
               5, 'MarkerFaceColor',[1 1 1]/4, 'MarkerEdgeColor',[1 1 1]/8)
            for k = 1:length(stripe_id_vec)
                stripe_id = stripe_id_vec(k);
                plot(spline_mat(:,k),1:size(stripe_mat,1),'Color',cm(1+(stripe_id-1)*increment,:),...
                        'LineWidth',1.5)
            end
            axis([0 size(stripe_mat,2) 0 size(stripe_mat,1)])
    %         text(10,25,[iIndex(round(t/60),2),...
    %         ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])
            ax = gca;
            ax.Visible = 'off';
            mkdir([flux_dynamics_path '/set_' num2str(i) '/'])
            saveas(centroid_fig,[flux_dynamics_path '/set_' num2str(i) '/ct_fluo_t' num2str(t) '.tif'],'tif');
        end
        temp_fluo_array(:,:,j) = frame_array;
        temp_spline_mat(:,:,j) = spline_mat;
        temp_disp_mat(:,:,j) = dispersion_mat;
    end    
    
    close all        
    %%%------------- Make time-dependent inference regions -------------%%%
    mean_fluo_mat = nanmean(temp_fluo_array,3);
    mean_center_mat = round(nanmean(temp_spline_mat,3));
    stripe_id_mat_full = NaN(size(stripe_mat,1),size(stripe_mat,2),length(t_vec));         
    for t = 1:length(t_vec)        
        center_mat = round(temp_spline_mat(:,:,t));
        for j = 1:length(stripe_id_vec)
            for k = 1:size(stripe_mat,1)
                stripe_id_mat_full(k,center_mat(k,j)-stripe_radius:...
                                center_mat(k,j)+stripe_radius,t) = stripe_id_vec(j);
            end  
            if j == 1 %&& j ~= length(stripe_id_vec)
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,max(1,center_mat(k,j) - 3*stripe_radius):...
                                center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                round(.5*(center_mat(k,j)+center_mat(k,j+1))),t) =...
                                stripe_id_vec(j) + 1/3;                            
                end                                                 
            elseif j == length(stripe_id_vec)
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                min(1024,center_mat(k,j) + 3*stripe_radius),t) =...
                                stripe_id_vec(j) + 1/3;                             
                end     
            else
                for k = 1:size(stripe_mat,1)
                    stripe_id_mat_full(k,round(.5*(center_mat(k,j)+center_mat(k,j-1))):...
                                    center_mat(k,j)-stripe_radius-1,t) = stripe_id_vec(j) - 1/3;
                    stripe_id_mat_full(k,min(1024,center_mat(k,j) + stripe_radius + 1):...
                                    round(.5*(center_mat(k,j)+center_mat(k,j+1))),t) =...
                                    stripe_id_vec(j) + 1/3;                            
                end
            end
        end
    end
    mode_stripe_id_mat = mode(stripe_id_mat_full,3);
    mean_stripe_id_mat = mean(stripe_id_mat_full,3);
    stripe_pos_struct(i).stripe_id_mat = stripe_id_mat_full;
     
    stripe_RGB = NaN(size(stripe_mat,1),size(stripe_mat,2),3);    
    index_vec = 1:size(stripe_mat,2);        
    % make RGB image
    for j = 1:length(stripe_id_vec)
        stripe_id = stripe_id_vec(j);
        stripe_color = cm(1+(stripe_id-1)*increment,:);
        for k = 1:3
            slice = stripe_RGB(:,:,k);
            slice(mode_stripe_id_mat==stripe_id-1/3) = stripe_color(k)/2;
            slice(mode_stripe_id_mat==stripe_id) = stripe_color(k)/1.5;
            slice(mode_stripe_id_mat==stripe_id+1/3) = stripe_color(k);
            slice(isnan(slice)) = .5;
            stripe_RGB(:,:,k) = slice;
        end
    end
    
    partition_fig = figure;  
    partition_fig.Visible = 'off';
    imshow(stripe_RGB)
    hold on
    for j = 1:length(stripe_id_vec)
        plot(mean_center_mat(:,j),1:size(stripe_mat,1),'Color','black','LineWidth',1.5)
    end
    title(['Set: ' num2str(i)])       
    saveas(partition_fig,[flux_dynamics_path '/inference_partitions_set' num2str(i) '.png'],'png') 
    
    %%% classify traces
    single_set_vec = [trace_struct.setID];
    set_indices = find(single_set_vec==i);
    
    for m = 1:length(set_indices)
        ind = set_indices(m);
        ParticleID = trace_struct(ind).ParticleID;
        xVec = xPos_vec_particle(particle_vec==ParticleID);
        yVec = yPos_vec_particle(particle_vec==ParticleID);
        t_trace = trace_struct(ind).time;
        f_trace = trace_struct(ind).fluo;
        t_trace = t_trace(~isnan(f_trace));
        tr_stripe_id_vec = zeros(1,length(t_trace))-1;
        for t = 1:length(t_trace)
            xp = xVec(t);
            yp = yVec(t);
            filter = round(t_trace(t)/60)==t_vec;
            if sum(filter) == 0
                continue
            end
            s_id = stripe_id_mat_full(yp,xp,filter);
            if ~isnan(s_id)
                dispersion_factor = nanmean(temp_disp_mat(:,stripe_id_vec==round(s_id),filter));
                if dispersion_factor <= max_disp
                    tr_stripe_id_vec(t) = s_id;
                else
                    tr_stripe_id_vec(t) = -1;
                end
            else
                tr_stripe_id_vec(t) = NaN;
            end
        end         
        fill_ind = find(-1==(tr_stripe_id_vec));
        id_vec = 1:length(tr_stripe_id_vec);
        bad_times = round(t_trace(fill_ind)/60);
        t_ref_vec = t_vec;
        t_ref_vec(ismember(t_vec,bad_times)) = Inf;
        if ~isempty(fill_ind) 
            % Assign early time points to nearest valid dynamic bin        
            for k = 1:length(fill_ind)
                [~, nn] = min(abs(t_ref_vec - round(t_trace(fill_ind(k))/60))); % find nearest valid point
                tr_stripe_id_vec(fill_ind(k)) = stripe_id_mat_full(yVec(fill_ind(k)),xVec(fill_ind(k)),nn);
            end            
        end        
        trace_struct(ind).stripe_id_inf = mode(tr_stripe_id_vec);
        trace_struct(ind).stripe_id_vec = tr_stripe_id_vec;
    end    
    
    %%% classify nuclei
    single_set_vec = [schnitz_struct.setID];
    set_indices = find(single_set_vec==i);
    
    for m = 1:length(set_indices)
        ind = set_indices(m);
        ncID = schnitz_struct(ind).ncID;
        xVec = round(xPos_vec_nc(nc_vec==ncID));
        yVec = round(yPos_vec_nc(nc_vec==ncID));
        t_trace = schnitz_struct(ind).time;                
        nc_stripe_id_vec = NaN(1,length(t_trace));
        for t = 1:length(t_trace)
            xp = xVec(t);
            yp = yVec(t);
            filter = round(t_trace(t)/60)==t_vec;
            if sum(filter) == 0
                continue
            end
            s_id = stripe_id_mat_full(yp,xp,filter);
            if ~isnan(s_id)
                dispersion_factor = nanmean(temp_disp_mat(:,stripe_id_vec==round(s_id),filter));
                if dispersion_factor <= max_disp
                    nc_stripe_id_vec(t) = s_id;
                else
                    nc_stripe_id_vec(t) = -1;
                end
            else
                nc_stripe_id_vec(t) = NaN;
            end
        end         
        fill_ind = find(-1==(nc_stripe_id_vec));
        id_vec = 1:length(nc_stripe_id_vec);
        bad_times = round(t_trace(fill_ind)/60);
        t_ref_vec = t_vec;
        t_ref_vec(ismember(t_vec,bad_times)) = Inf;
        if ~isempty(fill_ind) 
            % Assign early time points to nearest valid dynamic bin        
            for k = 1:length(fill_ind)
                [~, nn] = min(abs(t_ref_vec - round(t_trace(fill_ind(k))/60))); % find nearest valid point
                nc_stripe_id_vec(fill_ind(k)) = stripe_id_mat_full(yVec(fill_ind(k)),xVec(fill_ind(k)),nn);
            end            
        end                
        schnitz_struct(ind).stripe_id_inf = mode(nc_stripe_id_vec);
        schnitz_struct(ind).stripe_id_vec = nc_stripe_id_vec;
    end
    stripe_pos_struct(i).t_vec = t_vec;
    disp(['Completed ' num2str(i) ' of ' num2str(length(set_index))])    
end
save(stripe_save_name, 'stripe_pos_struct')
save(trace_name,'trace_struct')
save(nucleus_name,'schnitz_struct')