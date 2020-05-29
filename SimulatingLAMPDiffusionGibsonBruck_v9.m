%SimulatingLAMPDiffusionGibsonBruck_v5
%Update from v1 removing H ion diffusion within DNA clusters.
%Update from v2 removes H ions when they come into contact with DNA
%clusters instead of reflecting. Plotting of total H ions at all sensors.
%Prettied up figures.
%Update from v3 includes standard deviation of sensor output.
%Update from v4 to incorporate changes made to initialising script
%including realistic chip layouts and electrodes to manipulate DNA starting
%position.
%Update from v5 adding surface reactions for ISFET surface
%Update from v6 improving plotting
%Update from v7 making sensors work with all subvolume sizes. Added input
%checks
%Update from v8 moving chip building and inut checks to the initialisation code

addpath('LewisFunctions/');

clusters = clusters_init;

% Keep track of actual # of protons released
%tot_protons = 0;
nH = 0;

%Track total number of H ions at all sensors, for plotting over course of
%reaction, as well as number of H ions at each individual sensor
sensor_nH = zeros(N_x*N_y, length(t));
allsensor_nH = zeros(1, length(t));
plot_sensor_nH = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2));

%track standard deviation of sensor ouput over time
sensor_sd = zeros(1, length(t));

%Subvolume_nH is a 3D array designed to contain the number of H ions in
%each subvolume. Initialised to zero.
subvolume_nH = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));

%Create 3D array that contains reaction times for diffusion in each
%subvolume. Since starting with no H ions in array, all reaction times are
%infinity. Note these are absolute times instead of time to next reaction.
reactiontimes_diffusion = Inf(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));

%Create 2D array to track H ions bound to surface sites
surface_nH = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2));

%Create 2D array to store reaction times for surface binding
reactiontimes_surface = inf(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2));

%Create 3D array tracking the subvolumes within DNA clusters ( 0 = not in
%cluster, 1 = in cluster)
subvolume_incluster = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));

%Plotflag = 1 when a H ion has been added to the simulation. To avoid
%unnecessary plotting of figures
plotflag = 0;

%visualiseplot_flag = 0 sets plot visibility to 'off' to improve speed,
%particularly in remote desktop. visualiseplot_flag = 1 allows plots to be
%viewed as they are produced
visualiseplot_flag = 0;

%Create separate time variable to keep track of time during H ion diffusion
%simulation
t_H = 0;

for i = 2:length(t)
    
    num_of_active_sig = sum(i>sigmoid_start);
    cla(gca);
    
    % For every cluster...
    for j = 1:num_of_active_sig
        
        % If cluster is active then update propogation front
        prev_r_h = clusters(j).radius_h;
        prev_r_v = clusters(j).radius_v;
        
        clusters(j).radius_h = speed_of_front*(t(i)-t(sigmoid_start(j))); %BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, sensor_dim, N);
        clusters(j).radius_v = speed_of_front*(t(i)-t(sigmoid_start(j))); %BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, solution_height, M);
        
        % Compute position of randomly distributed protons around propogation front
        numtorelease = round(burst*(volume_of_sphere(clusters(j).radius_h)  - volume_of_sphere(prev_r_h)));
        
        c = 2*rand(numtorelease,1)-1;
        lon=2*pi*rand(numtorelease,1);
        lat=acos(c);
        a=cos(lon).*sin(lat);
        b=sin(lon).*sin(lat);
        
        additions = [clusters(j).radius_h*a'+clusters(j).centre_x;clusters(j).radius_h*b'+clusters(j).centre_y;clusters(j).radius_v*c'+clusters(j).centre_z]';
        
        % Remove the computed proton positions outside the array
        additions(additions(:,1)>sensor_xsize,:) = [];
        additions(additions(:,1)<0,:) = [];
        additions(additions(:,2)>sensor_ysize,:) = [];
        additions(additions(:,2)<0,:) = [];
        additions(additions(:,3)>solution_height,:) = [];
        additions(additions(:,3)<0,:) = [];
        
        for k = 1:num_of_mol
            if(k~=j)
                additions( (additions(:,1) - clusters(k).centre_x).^2 ...
                    + (additions(:,2) - clusters(k).centre_y).^2 ...
                    + (additions(:,3) - clusters(k).centre_z).^2 ...
                    <= clusters(k).radius_h^2, : ) = [];
            end
        end
        
        %Calculate distance from each subvolume to cluster surface
        subvolumetoclusterradius_mag = ((subvolume_mesh3x - clusters(j).centre_x).^2 +(subvolume_mesh3y - clusters(j).centre_y).^2 +(subvolume_mesh3z - clusters(j).centre_z).^2).^0.5 - clusters(j).radius_h;        
        
        %Find index of subvolumetoclusterradius values that are negative.
        %Note meshgrid swaps x and y dimensions
        subvolumeincluster_linear = find(subvolumetoclusterradius_mag < 0);
        [subvolumeincluster_indexy, subvolumeincluster_indexx, subvolumeincluster_indexz] = ind2sub(size(subvolumetoclusterradius_mag), subvolumeincluster_linear);
        
        if size(subvolumeincluster_linear, 1) > 0
            
            %Remove H ions from subvolumes inside cluster
            for l = 1:length(subvolumeincluster_linear)

                subvolume_incluster(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) = 1;

                %Remove H ions within expanded cluster and set reaction
                %time to inf
                if subvolume_nH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) > 0  

                    %Set reactiontime = inf and nH = 0 if subvolume within
                    %cluster
                    reactiontimes_diffusion(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) = inf;
                    reactiontimes_surface(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l)) = inf;
                    subvolume_nH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) = 0;
                end

            end            
        end
        
        %Find the subvolume which the added H ions belong to and add them
        %to subvolume_nH
        if size(additions, 1) > 0
            
            if plotflag == 0
                plotflag = 1;
                t_H = t(i)
            end
                        
            for k = 1:size(additions, 1)
                
                %Calculate distance between each subvolume and the addition.
                %Set values for subvolumes within cluster to inf, so that 
                %min function finds nearest subvolume outside cluster
                subvolumetoaddition_mag = ((subvolume_mesh3x - additions(k, 1)).^2 +(subvolume_mesh3y - additions(k, 2)).^2 +(subvolume_mesh3z - additions(k, 3)).^2).^0.5;
                if ~isempty(subvolumeincluster_linear)
                    for l = 1:length(subvolumeincluster_linear)
                        %Note subvolumetoaddition_mag has dimensions y, x,
                        %z as meshgrid swaps x and y round
                        subvolumetoaddition_mag(subvolumeincluster_indexy(l), subvolumeincluster_indexx(l), subvolumeincluster_indexz(l)) = inf;
                    end
                end
                
                %Find nearest subvolume. Note meshgrid swaps x and y
                %dimensions
                [~, additionsubvolume_linear] = min(subvolumetoaddition_mag(:));
                [additionsubvolume_indexy, additionsubvolume_indexx, additionsubvolume_indexz] = ind2sub(size(subvolumetoaddition_mag), additionsubvolume_linear);

                %Add H ion to subvolume_nH
                subvolume_nH(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz) = subvolume_nH(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz) + 1;
                nH = nH + 1;
                
                %Compute new diffusion reaction time for subvolume H ion has been
                %added to (See equation 3.9 in
                %https://core.ac.uk/download/pdf/1568321.pdf).
                reactiontimes_diffusion(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz) = ReactionTimeRecalculation_nHChange( additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz, subvolume_nH, reactiontimes_diffusion(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz), k_diffusion, t_H);
                
                %Compute new surface reaction time if H ion added to z=0
                %plane
                if additionsubvolume_indexz == 1
                    reactiontimes_surface(additionsubvolume_indexx, additionsubvolume_indexy) = ReactionTimeRecalculation_nHChange( additionsubvolume_indexx, additionsubvolume_indexy, 1, subvolume_nH, reactiontimes_surface(additionsubvolume_indexx, additionsubvolume_indexy), k_surface, t_H);
                end
            end
        end
    end
    
    %Simulate H ion diffusion if they exist within solution and not in z=0
    %plane, where H ions are assumed to be fixed to ISFET
    %Carry out diffusion of H ions until next time step for LAMP
    %reaction takes place    
    while i ~= size(t,2) && t_H < t(i+1) && nH ~= 0
        
        %Find smallest diffusion event time and surface reaction event time
        %and determine which is smallest
        [diffusion_time, diffusion_linear] = min(reactiontimes_diffusion(:));
        [surface_time, surface_linear] = min(reactiontimes_surface(:));
        
        if diffusion_time < surface_time
            
            event_time = diffusion_time;
            [event_indexx, event_indexy, event_indexz] = ind2sub(size(reactiontimes_diffusion), diffusion_linear);
        
            %Draw random variable event_type
            %between 1-6, corresponding to +x, -x, +y, -y, +z, -z, diffusion respectively
            %. Possible as all directions of diffusion are equally likely.
            event_type = randi(6);
        
        else
            event_time = surface_time;
            [event_indexx, event_indexy] = ind2sub(size(reactiontimes_surface), surface_linear);
            event_type = 7;
        end
        
        %Only carry out event if it falls before the next step in the LAMP
        %reaction
        if event_time < t(i+1)          
            
            %Update time
            t_H = event_time
            
            if event_type <=6
                
                %No diffusion beyond boundaries of the reaction chamber. H
                %ion stays in same subvolume, and reactiontime is
                %recalculated
                if event_indexx == size(subvolume_nH, 1) && event_type == 1 || event_indexx == 1 && event_type == 2 || event_indexy == size(subvolume_nH, 2) && event_type == 3 || event_indexy == 1 && event_type == 4 || event_indexz == size(subvolume_nH, 3) && event_type == 5 || event_indexz == 1 && event_type == 6
                    disp('Hit wall')
                    reactiontimes_diffusion(event_indexx, event_indexy, event_indexz) = ReactionTimeRecalculation_EventOccurred( event_indexx, event_indexy, event_indexz, subvolume_nH, k_diffusion, t_H);

                %No diffusion if moving into subvolume within DNA cluster.
                %H ion stays in same subvolume, and reactiontime is
                %recalculated
                elseif event_type == 1 && subvolume_incluster(event_indexx+1, event_indexy, event_indexz) == 1 || event_type == 2 && subvolume_incluster(event_indexx-1, event_indexy, event_indexz) == 1 || event_type == 3 && subvolume_incluster(event_indexx, event_indexy+1, event_indexz) == 1 || event_type == 4 && subvolume_incluster(event_indexx, event_indexy-1, event_indexz) == 1 || event_type == 5 && subvolume_incluster(event_indexx, event_indexy, event_indexz+1) == 1 || event_type == 6 && subvolume_incluster(event_indexx, event_indexy, event_indexz-1) == 1 
                    disp('Hit DNA cluster')
                    
                    %Remove H ion from simulation if
                    %Hclusterinteraction_flag = 1
                    if Hclusterinteraction_flag == 1
                       subvolume_nH(event_indexx, event_indexy, event_indexz) = subvolume_nH(event_indexx, event_indexy, event_indexz) - 1;
                       
                       %If H ion removed from z=0 plane reclaculate surface
                       %reaction time
                       if event_indexz == 1
                            reactiontimes_surface(event_indexx, event_indexy) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy, 1, subvolume_nH, reactiontimes_surface(event_indexx, event_indexy), k_surface, t_H);
                       end
                    end
                    
                    reactiontimes_diffusion(event_indexx, event_indexy, event_indexz) = ReactionTimeRecalculation_EventOccurred( event_indexx, event_indexy, event_indexz, subvolume_nH, k_diffusion, t_H);
                    
                else
                    
                    %Remove H ion from subvolume, and recalculate reaction
                    %time
                    subvolume_nH(event_indexx, event_indexy, event_indexz) = subvolume_nH(event_indexx, event_indexy, event_indexz) - 1;
                    reactiontimes_diffusion(event_indexx, event_indexy, event_indexz) = ReactionTimeRecalculation_EventOccurred( event_indexx, event_indexy, event_indexz, subvolume_nH, k_diffusion, t_H);
                    
                    %Recalculate surface reaction time if H ion was
                    %originally in z=0 plane
                    if event_indexz == 1
                        reactiontimes_surface(event_indexx, event_indexy) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy, 1, subvolume_nH, reactiontimes_surface(event_indexx, event_indexy), k_surface, t_H);
                    end
                    
                    %Carry out diffusion event and recalculate reaction
                    %time
                    if event_type == 1
                        %Diffusion +x
                        disp('Diffusion +x')

                        subvolume_nH(event_indexx+1, event_indexy, event_indexz) = subvolume_nH(event_indexx+1, event_indexy, event_indexz) + 1;
                        reactiontimes_diffusion(event_indexx+1, event_indexy, event_indexz) = ReactionTimeRecalculation_nHChange( event_indexx+1, event_indexy, event_indexz, subvolume_nH, reactiontimes_diffusion(event_indexx+1, event_indexy, event_indexz), k_diffusion, t_H);
                        
                        %Recalculate surface reaction time if H ion moving
                        %into z=0 plane
                        if event_indexz == 1
                            reactiontimes_surface(event_indexx+1, event_indexy) = ReactionTimeRecalculation_nHChange( event_indexx+1, event_indexy, 1, subvolume_nH, reactiontimes_surface(event_indexx+1, event_indexy), k_surface, t_H);
                        end
                        
                    elseif event_type == 2
                        %Diffusion -x
                        disp('Diffusion - x')         

                        subvolume_nH(event_indexx-1, event_indexy, event_indexz) = subvolume_nH(event_indexx-1, event_indexy, event_indexz) + 1;
                        reactiontimes_diffusion(event_indexx-1, event_indexy, event_indexz) = ReactionTimeRecalculation_nHChange( event_indexx-1, event_indexy, event_indexz, subvolume_nH, reactiontimes_diffusion(event_indexx-1, event_indexy, event_indexz), k_diffusion, t_H);
                        
                        %Recalculate surface reaction time if H ion moving
                        %into z=0 plane
                        if event_indexz == 1
                            reactiontimes_surface(event_indexx-1, event_indexy) = ReactionTimeRecalculation_nHChange( event_indexx-1, event_indexy, 1, subvolume_nH, reactiontimes_surface(event_indexx-1, event_indexy), k_surface, t_H);
                        end
                        
                    elseif event_type == 3
                        %Diffusion +y
                        disp('Diffusion +y')

                        subvolume_nH(event_indexx, event_indexy+1, event_indexz) = subvolume_nH(event_indexx, event_indexy+1, event_indexz) + 1;
                        reactiontimes_diffusion(event_indexx, event_indexy+1, event_indexz) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy+1, event_indexz, subvolume_nH, reactiontimes_diffusion(event_indexx, event_indexy+1, event_indexz), k_diffusion, t_H);
                        
                        %Recalculate surface reaction time if H ion moving
                        %into z=0 plane
                        if event_indexz == 1
                            reactiontimes_surface(event_indexx, event_indexy+1) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy+1, 1, subvolume_nH, reactiontimes_surface(event_indexx, event_indexy+1), k_surface, t_H);
                        end
                        
                    elseif event_type == 4
                        %Diffusion -y
                        disp('Diffusion -y')

                        subvolume_nH(event_indexx, event_indexy-1, event_indexz) = subvolume_nH(event_indexx, event_indexy-1, event_indexz) + 1;
                        reactiontimes_diffusion(event_indexx, event_indexy-1, event_indexz) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy-1, event_indexz, subvolume_nH, reactiontimes_diffusion(event_indexx, event_indexy-1, event_indexz), k_diffusion, t_H);
                        
                        %Recalculate surface reaction time if H ion moving
                        %into z=0 plane
                        if event_indexz == 1
                            reactiontimes_surface(event_indexx, event_indexy-1) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy-1, 1, subvolume_nH, reactiontimes_surface(event_indexx, event_indexy-1), k_surface, t_H);
                        end
                        
                    elseif event_type == 5
                        %Diffusion +z
                        disp('Diffusion +z')

                        subvolume_nH(event_indexx, event_indexy, event_indexz+1) = subvolume_nH(event_indexx, event_indexy, event_indexz+1) + 1;
                        reactiontimes_diffusion(event_indexx, event_indexy, event_indexz+1) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy, event_indexz+1, subvolume_nH, reactiontimes_diffusion(event_indexx, event_indexy, event_indexz+1), k_diffusion, t_H);
                        
                    elseif event_type == 6
                        %Diffusion -z
                        disp('Diffusion -z')

                        subvolume_nH(event_indexx, event_indexy, event_indexz-1) = subvolume_nH(event_indexx, event_indexy, event_indexz-1) + 1;
                        reactiontimes_diffusion(event_indexx, event_indexy, event_indexz-1) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy, event_indexz-1, subvolume_nH, reactiontimes_diffusion(event_indexx, event_indexy, event_indexz-1), k_diffusion, t_H);
                        
                        %Recalculate surface reaction time if H ion moving
                        %into z=0 plane
                        if event_indexz == 2
                            reactiontimes_surface(event_indexx, event_indexy) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy, 1, subvolume_nH, reactiontimes_surface(event_indexx, event_indexy), k_surface, t_H);
                        end
                        
                    end     
                end

            elseif event_type == 7
                %Surface Reaction
                disp('Surface reaction')
                
                %Remove H ion from z=0 plane of subvolume_nH and
                %recalculate both the diffusion and surface reaction times
                subvolume_nH(event_indexx, event_indexy, 1) = subvolume_nH(event_indexx, event_indexy, 1) - 1;
                reactiontimes_diffusion(event_indexx, event_indexy, 1) = ReactionTimeRecalculation_nHChange( event_indexx, event_indexy, 1, subvolume_nH, reactiontimes_diffusion(event_indexx, event_indexy, 1), k_diffusion, t_H);
                reactiontimes_surface(event_indexx, event_indexy) = ReactionTimeRecalculation_EventOccurred( event_indexx, event_indexy, 1, subvolume_nH, k_surface, t_H);
            
                %Add H ion to surface site
                surface_nH(event_indexx, event_indexy) = surface_nH(event_indexx, event_indexy) + 1;
                
            end
        
        %If next reaction time is after next LAMP reaction step, set the
        %simulation time to the next LAMP time step
        else
            t_H = t(i+1)
        end
        
        if t_H == -inf
            return
        end
        
    end
    
    %Plot figures if plotflag = 1 (H ions exist in simulation)
    if plotflag == 1

        %Extract data for plotting
        plot_nH_linear = find(subvolume_nH > 0);
        [plot_nH_indexx, plot_nH_indexy, plot_nH_indexz] = ind2sub( size(subvolume_nH), plot_nH_linear);
        plot_x_nH = subvolume_xyz{1,1}(plot_nH_indexx);
        plot_y_nH = subvolume_xyz{1,2}(plot_nH_indexy);
        plot_z_nH = subvolume_xyz{1,3}(plot_nH_indexz);
        plot_nH = subvolume_nH(plot_nH_linear);
        
        plot_incluster_linear = find(subvolume_incluster > 0);
        [plot_incluster_indexx, plot_incluster_indexy, plot_incluster_indexz] = ind2sub( size(subvolume_incluster), plot_incluster_linear);
        plot_x_incluster = subvolume_xyz{1,1}(plot_incluster_indexx);
        plot_y_incluster = subvolume_xyz{1,2}(plot_incluster_indexy);
        plot_z_incluster = subvolume_xyz{1,3}(plot_incluster_indexz);
        plot_incluster = subvolume_incluster(plot_incluster_linear);
        
        clf
        
        %Plot 3D scatter plot of proton locations and save
        figure(1)
        scatter3(plot_x_nH, plot_y_nH, plot_z_nH, plot_nH)
        if visualiseplot_flag == 0
            set(gcf, 'Visible', 'off')
        end
        axis([0 sensor_xsize 0 sensor_ysize 0 solution_height])
        title(sprintf('Proton Locations (t = %d)', t_H))
        xlabel('x'); ylabel('y'); zlabel('z');
        hold on
        h = surface([0, 0; sensor_xsize, sensor_xsize], [0, sensor_ysize; 0, sensor_ysize], [0, 0; 0, 0], chip_layout,'facecolor', 'texturemap', 'edgecolor', 'none');
        alpha(h, 0.5);
        hold off
        if visualiseplot_flag == 1
            shg     
        end
        filename1 = sprintf('%s\\ProtonLocations_t%d.fig', save_D, int16(t(i)));
        savefig(filename1)
       
        %Plot 3D scatter plot of DNA clusters and save
        figure(2)
        scatter3(plot_x_incluster, plot_y_incluster, plot_z_incluster, plot_incluster)
        if visualiseplot_flag == 0
            set(gcf, 'Visible', 'off')
        end
        axis([0 sensor_xsize 0 sensor_ysize 0 solution_height])
        title(sprintf('DNA Clusters (t = %d)', t_H))         
        xlabel('x'); ylabel('y'); zlabel('z');
        hold on
        h = surface([0, 0; sensor_xsize, sensor_xsize], [0, sensor_ysize; 0, sensor_ysize], [0, 0; 0, 0], chip_layout,'facecolor', 'texturemap', 'edgecolor', 'none');
        alpha(h, 0.5);
        hold off
        if visualiseplot_flag == 1
            shg        
        end
        filename2 = sprintf('%s\\Clusters_t%d.fig', save_D, int16(t(i)));
        savefig(filename2)
        
        %Sum up H ions at z=0 plane into each sensor. j represents the
        %sensor (row of sensor_location).
        for j = 1:size(sensor_location,1)
            sensor_nH_temp = 0;
            for k = sensor_location{j, 1}
                sensor_nH_temp = sensor_nH_temp + surface_nH(k);
            end
            sensor_nH(j, i) = sensor_nH_temp;
            
            %Save sensor_nH values for surface plotting
            for k = sensor_location{j, 1}
                plot_sensor_nH(k) = sensor_nH(j, i);
            end
        end
        
        figure(3)
        surface(subvolume_mesh2x, subvolume_mesh2y, plot_sensor_nH.')
        if visualiseplot_flag == 0
            set(gcf, 'Visible', 'off')
        end
        title(sprintf('Sensor #H (t = %d)', t_H))
        xlabel('x'); ylabel('y'); zlabel('# H+');
        [z_max, ~] = max(plot_sensor_nH(:));
        if z_max == 0
            z_max = 1;
        end
        axis([0 sensor_xsize 0 sensor_ysize 0 z_max]);
        if visualiseplot_flag == 1
            shg        
        end
        filename3 = sprintf('%s\\SensorH_t%d.fig', save_D, int16(t(i)));
        savefig(filename3)

        %Compute standard deviation, save value, plot line graph
        %{
        sensor_sd(1, i) = std(current_sensor_nH, 0, 'all');
        figure(5)
        plot(t, sensor_sd)
        title(sprintf('Sensor Standard Deviation (t = %d)', t_H))
        shg
        %}
                
        %Save graph if end of simulation
        %if i == length(t)-1
            %filename5 = sprintf('%s\\SensorSD_t%d.fig', save_D, int16(t(i)));
            %savefig(filename5)
        %end
        
        %Sum nH cumulatively for all sensors, averaging over number of
        %sensors
        allsensor_nH(1, i) = sum(sensor_nH(:, i))/(N_x*N_y);
        
        %Sum H ions over all sensors, save value, plot line graph
        figure(4)
        plot(t, allsensor_nH)
        if visualiseplot_flag == 0
            set(gcf, 'Visible', 'off')
        end
        title(sprintf('Sensor #H average (t = %d)', t_H))
        xlabel('Time (s)'); ylabel('Average # H+');
        if visualiseplot_flag == 1
            shg        
        end
        
        %Save graph if end of simulation
        if i == length(t)-1
            filename4 = sprintf('%s\\SensorAverage_t%d.fig', save_D, int16(t(i)));
            savefig(filename4)
        end
        
    end
    
end
