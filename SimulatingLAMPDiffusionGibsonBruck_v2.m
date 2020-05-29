%SimulatingLAMPDiffusionGibsonBruck_v3
%Update from v1 removing H ion diffusion within DNA clusters

clusters = clusters_init;

% Create figure
f = figure;

Z = zeros(N, N);
track = zeros(length(t), N, N);

persisted_protons = 0;

% Keep track of actual # of protons released
tot_protons = 0;

s = surf(zeros(N,N));
xlim([1 N]);
ylim([1 N]);
zlim([1 N]);
hold on;

set(gca,'FontSize',16);
mystr = {'Release of Protons in 3D Space',['t = 0 minutes']};
    title(mystr, 'FontSize', 20);

xlabel('x', 'FontSize', 16)
ylabel('y', 'FontSize', 16)
zlabel('z', 'FontSize', 16)

%Create 3D array compartmentalising reaction chamber and store the number
%of H ions in each subvolume
subvolume_xyzsize = [1, 1, 1];
subvolume_xyz = {subvolume_xyzsize(1)/2:subvolume_xyzsize(1):N-subvolume_xyzsize(1)/2, subvolume_xyzsize(2)/2:subvolume_xyzsize(2):N-subvolume_xyzsize(2)/2, 0:subvolume_xyzsize(3):M};
[subvolume_meshx, subvolume_meshy, subvolume_meshz] = meshgrid(subvolume_xyz{1,1}, subvolume_xyz{1,2}, subvolume_xyz{1,3});
subvolume_nH = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));
[sensor_meshx, sensor_meshy] = meshgrid(0.5:N-0.5, 0.5:N-0.5);

subvolume_n = size(subvolume_xyz{1,1},2)*size(subvolume_xyz{1,2},2)*size(subvolume_xyz{1,3},2);
xstep = 1/subvolume_xyzsize(1);
ystep = 1/subvolume_xyzsize(2);
plot_meshx = reshape(subvolume_meshx, [1, subvolume_n]);
plot_meshy = reshape(subvolume_meshx, [1, subvolume_n]); 
plot_meshz = reshape(subvolume_meshx, [1, subvolume_n]); 

nH = 0; 
nH_incluster = 0;
plotflag = 0;

%Create 3D array tracking the subvolumes within DNA clusters ( 0 = not in
%cluster, 1 = in cluster)
subvolume_incluster = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));

%Create 3D array to track removed and added H ions
removedH = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));
addedH = zeros(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2));

%Create 3D array that contains reaction times for diffusion in each
%subvolume. Since starting with no H ions in array, all reaction times are
%infinity. Note these are absolute times instead of time to next reaction.
%-1 in z dimension to allow for no calculation of reaction times for z=0
%plane, where H ions associated with ISFET sensors
reactiontimes = Inf(size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2)-1);

%Create separate time variable to keep track of time during H ion diffusion
%simulation
t_H = 0;

%Effective 'rate constant' for H ion diffusion, k_diffusion.
%Value from http://omh.umeche.maine.edu/pdfs/JChemPhys_135_124505.01pdf.pdf
%for 50 degrees C. Units m^2/s.
diffusivity_H = 1E-8;
k_diffusion = diffusivity_H/(0.03125E-3)^2;

for i = 2:length(t)
    
    num_of_active_sig = sum(i>sigmoid_start);
    cla(gca);
    
    % For every cluster...
    for j = 1:num_of_active_sig
        % If cluster is active then update propogation front
        prev_r_h = clusters(j).radius_h;
        prev_r_v = clusters(j).radius_v;
        
        clusters(j).radius_h = BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, sensor_dim, N);
        clusters(j).radius_v = BuMeters(speed_of_front*(t(i)-t(sigmoid_start(j))), 1, solution_height, M);
        
        % Compute position of randomly distributed protons around propogation front
        numtorelease = round(burst*(volume_of_sphere(clusters(j).radius_h)  - volume_of_sphere(prev_r_h)));
        
        c = 2*rand(numtorelease,1)-1;
        lon=2*pi*rand(numtorelease,1);
        lat=acos(c);
        a=cos(lon).*sin(lat);
        b=sin(lon).*sin(lat);
        
        additions = [clusters(j).radius_h*a'+clusters(j).centre_x;clusters(j).radius_h*b'+clusters(j).centre_y;clusters(j).radius_v*c'+clusters(j).centre_z]';
        
        % Remove the computed proton positions outside the array
        additions(additions(:,1)>N+1,:) = [];
        additions(additions(:,1)<1,:) = [];
        additions(additions(:,2)>N+1,:) = [];
        additions(additions(:,2)<1,:) = [];
        additions(additions(:,3)>M+1,:) = [];
        additions(additions(:,3)<1,:) = [];
        
        for k = 1:num_of_mol
            if(k~=j)
                additions( (additions(:,1) - clusters(k).centre_x).^2 ...
                    + (additions(:,2) - clusters(k).centre_y).^2 ...
                    + (additions(:,3) - clusters(k).centre_z).^2 ...
                    <= clusters(k).radius_h^2, : ) = [];
            end
        end
        
        %Calculate distance from each subvolume to cluster surface
        subvolumetoclusterradius_mag = sqrt((subvolume_meshx - clusters(j).centre_x).^2 +(subvolume_meshy - clusters(j).centre_y).^2 +(subvolume_meshz - clusters(j).centre_z).^2) - clusters(j).radius_h;        
        
        %Find index of subvolumetoclusterradius values that are negative
        subvolumeincluster_linear = find(subvolumetoclusterradius_mag < 0);
        [subvolumeincluster_indexx, subvolumeincluster_indexy, subvolumeincluster_indexz] = ind2sub([size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2)], subvolumeincluster_linear);
        
        if size(subvolumeincluster_linear, 1) > 0
            
            %Remove H ions from subvolumes inside cluster
            for l = 1:size(subvolumeincluster_linear, 1)
                
                %Set incluster flag to 1 for subvolumes in cluster, apart
                %from in z=0 plane
                if subvolumeincluster_indexz(l) ~= 1
                    
                    subvolume_incluster(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) = 1;
                    
                    %Remove H ions within expanded cluster and set reaction
                    %time to inf
                    if subvolume_nH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) > 0  

                        %Set reactiontime = inf and nH = 0 if subvolume within
                        %cluster, except for subvolumes at z=0, to avoid wiping
                        %ISFET accumulation
                        reactiontimes(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)-1) = inf;
                        removedH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) = removedH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) + subvolume_nH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l));
                        subvolume_nH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) = 0;
                        nH = nH - subvolume_nH(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l))
                        disp('H ion removed')
                    end
                end
            end            
        end
        
        %Find the subvolume which the added H ions belong to and add them
        %to subvolume_nH
        if size(additions, 1) > 0
            
            disp('H ions added.')
            if plotflag == 0
                plotflag = 1;
                t_H = t(i)
            end
                        
            for k = 1:size(additions, 1)
                
                %Calculate distance between each subvolume and the addition.
                %Set values for subvolumes within cluster to inf, so that 
                %min function finds nearest subvolume outside cluster
                subvolumetoaddition_mag = sqrt((subvolume_meshx - additions(k, 1)).^2 +(subvolume_meshy - additions(k, 2)).^2 +(subvolume_meshz - additions(k, 3)).^2);
                if size(subvolumeincluster_linear, 1) > 0
                    for l = 1:size(subvolumeincluster_linear, 2)
                        subvolumetoaddition_mag(subvolumeincluster_indexx(l), subvolumeincluster_indexy(l), subvolumeincluster_indexz(l)) = inf;
                    end
                end
                
                %Find nearest subvolume
                [~, additionsubvolume_linear] = min(subvolumetoaddition_mag(:));
                [additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz] = ind2sub([size(subvolume_xyz{1,1},2), size(subvolume_xyz{1,2},2), size(subvolume_xyz{1,3},2)], additionsubvolume_linear);
                %disp(additionsubvolume_linear)
                %disp(additionsubvolume_indexx)
                %disp(subvolumetoclusterradius_mag)
                %disp(subvolumetoaddition_mag)
                %Add H ion to subvolume_nH
                subvolume_nH(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz) = subvolume_nH(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz) + 1;
                addedH(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz) = addedH(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz) + 1;
                nH = nH + 1
                
                %Compute new reaction time for subvolume H ion has been
                %added to (See equation 3.9 in
                %https://core.ac.uk/download/pdf/1568321.pdf).
                if additionsubvolume_indexz ~= 1
                    reactiontimes(additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz-1) = DiffusionReactionTimeRecalculation_In( additionsubvolume_indexx, additionsubvolume_indexy, additionsubvolume_indexz-1, subvolume_nH, reactiontimes, k_diffusion, t_H);
                end
                
            end
            
        end
    end
    
    %Simulate H ion diffusion if they exist within solution and not in z=0
    %plane, where H ions are assumed to be fixed to ISFET

    %Carry out diffusion of H ions until next time step for LAMP
    %reaction takes place
    
    while t_H < t(i+1) && nH ~= 0 && i ~= size(t,2)
        
        %Determine location of next diffusion event, with smallest reaction
        %time
        [event_time, event_linear] = min(reactiontimes(:));
        [event_indexx, event_indexy, event_indexz] = ind2sub(size(reactiontimes), event_linear);
        
        if event_time < t_H
            disp('Error. Reaction time earlier than current t_H')
        end
        %Carry out diffusion event, drawing random variable event_type
        %between 1-6, corresponding to +x, -x, +y, -y, +z, -z, respectively
        %. Possible as all directions of diffusion are equally likely.
        event_type = randi(6);
        
        %Note that reactiontimes is one element smaller than subvolume_nH
        %to skip reaction time calculation for H ions in the z=0 plane.
        %Therefore to access the corresponding subvolume in subvolume_nH
        %using the variable eventz, eventz+1 is the correct index.
        if event_time < t(i+1)          
            
            %Update time
            t_H = event_time
            
                if event_type <=6
                
                %No diffusion beyond boundaries of the reaction chamber. H
                %ion stays in same subvolume, and reactiontime is
                %recalculated
                if event_indexx == size(subvolume_nH, 1) && event_type == 1 || event_indexx == 1 && event_type == 2 || event_indexy == size(subvolume_nH, 2) && event_type == 3 || event_indexy == 1 && event_type == 4 || event_indexz == size(subvolume_nH, 3)-1 && event_type == 5 
                    disp('Hit wall')
                    reactiontimes(event_indexx, event_indexy, event_indexz) = DiffusionReactionTimeRecalculation_Out( event_indexx, event_indexy, event_indexz, subvolume_nH, k_diffusion, t_H);

                %No diffusion if moving into subvolume within DNA cluster.
                %H ion stays in same subvolume, and reactiontime is
                %recalculated
                elseif event_type == 1 && subvolume_incluster(event_indexx+1, event_indexy, event_indexz+1) == 1 || event_type == 2 && subvolume_incluster(event_indexx-1, event_indexy, event_indexz+1) == 1 || event_type == 3 && subvolume_incluster(event_indexx, event_indexy+1, event_indexz+1) == 1 || event_type == 4 && subvolume_incluster(event_indexx, event_indexy-1, event_indexz+1) == 1 || event_type == 5 && subvolume_incluster(event_indexx, event_indexy, event_indexz+2) == 1 || event_type == 6 && subvolume_incluster(event_indexx, event_indexy, event_indexz) == 1 
                    disp('Hit DNA cluster')                    
                    reactiontimes(event_indexx, event_indexy, event_indexz) = DiffusionReactionTimeRecalculation_Out( event_indexx, event_indexy, event_indexz, subvolume_nH, k_diffusion, t_H);
                    
                else
                    
                    %Remove H ion from subvolume, and recalculate reaction
                    %time
                    subvolume_nH(event_indexx, event_indexy, event_indexz+1) = subvolume_nH(event_indexx, event_indexy, event_indexz+1) - 1;
                    reactiontimes(event_indexx, event_indexy, event_indexz) = DiffusionReactionTimeRecalculation_Out( event_indexx, event_indexy, event_indexz, subvolume_nH, k_diffusion, t_H);
                    
                    %Carry out diffusion event and recalculate reaction
                    %time
                    if event_type == 1
                        %Diffusion +x
                        disp('Diffusion +x')

                        subvolume_nH(event_indexx+1, event_indexy, event_indexz+1) = subvolume_nH(event_indexx+1, event_indexy, event_indexz+1) + 1;
                        reactiontimes(event_indexx+1, event_indexy, event_indexz) = DiffusionReactionTimeRecalculation_In( event_indexx+1, event_indexy, event_indexz, subvolume_nH, reactiontimes, k_diffusion, t_H);
                        
                    elseif event_type == 2
                        %Diffusion -x
                        disp('Diffusion - x')         

                        subvolume_nH(event_indexx-1, event_indexy, event_indexz+1) = subvolume_nH(event_indexx-1, event_indexy, event_indexz+1) + 1;
                        reactiontimes(event_indexx-1, event_indexy, event_indexz) = DiffusionReactionTimeRecalculation_In( event_indexx-1, event_indexy, event_indexz, subvolume_nH, reactiontimes, k_diffusion, t_H);
                        
                    elseif event_type == 3
                        %Diffusion +y
                        disp('Diffusion +y')

                        subvolume_nH(event_indexx, event_indexy+1, event_indexz+1) = subvolume_nH(event_indexx, event_indexy+1, event_indexz+1) + 1;
                        reactiontimes(event_indexx, event_indexy+1, event_indexz) = DiffusionReactionTimeRecalculation_In( event_indexx, event_indexy+1, event_indexz, subvolume_nH, reactiontimes, k_diffusion, t_H);
                        
                    elseif event_type == 4
                        %Diffusion -y
                        disp('Diffusion -y')

                        subvolume_nH(event_indexx, event_indexy-1, event_indexz+1) = subvolume_nH(event_indexx, event_indexy-1, event_indexz+1) + 1;
                        reactiontimes(event_indexx, event_indexy-1, event_indexz) = DiffusionReactionTimeRecalculation_In( event_indexx, event_indexy-1, event_indexz, subvolume_nH, reactiontimes, k_diffusion, t_H);
                        
                    elseif event_type == 5
                        %Diffusion +z
                        disp('Diffusion +z')

                        subvolume_nH(event_indexx, event_indexy, event_indexz+2) = subvolume_nH(event_indexx, event_indexy, event_indexz+2) + 1;
                        reactiontimes(event_indexx, event_indexy, event_indexz+1) = DiffusionReactionTimeRecalculation_In( event_indexx, event_indexy, event_indexz+1, subvolume_nH, reactiontimes, k_diffusion, t_H);
                        
                    elseif event_type == 6
                        %Diffusion -z
                        disp('Diffusion -z')

                        subvolume_nH(event_indexx, event_indexy, event_indexz) = subvolume_nH(event_indexx, event_indexy, event_indexz) + 1;
                        
                        if event_indexz == 1
                            %Do not calculate reaction time if H ion is
                            %moving into z = 0 plane
                            nH = nH-1;                           
                        else
                            reactiontimes(event_indexx, event_indexy, event_indexz-1) = DiffusionReactionTimeRecalculation_In( event_indexx, event_indexy, event_indexz-1, subvolume_nH, reactiontimes, k_diffusion, t_H);
                        end
                        
                    end
                end

            else
                %Reaction
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
        
    plot_x_nH = [];
    plot_y_nH = [];
    plot_z_nH = [];
    plot_nH = [];
    
    plot_incluster = [];
    plot_x_incluster = [];
    plot_y_incluster = [];
    plot_z_incluster = [];
   
    plot_removedH = [];
    plot_x_removedH = [];
    plot_y_removedH = [];
    plot_z_removedH = [];
   
    for k = 1:size(subvolume_nH, 1)
        for l = 1:size(subvolume_nH, 2)
            for m = 1:size(subvolume_nH, 3)
                if subvolume_nH(k, l, m) > 0 && subvolume_incluster(k, l, m) == 1
                    disp('Error. H ion in DNA cluster.')
                    nH_incluster = nH_incluster + subvolume_nH(k, l, m);
                    %return
                end
                if subvolume_nH(k, l, m) > 0
                    plot_x_nH = [plot_x_nH, subvolume_xyz{1,1}(k)];
                    plot_y_nH = [plot_y_nH, subvolume_xyz{1,2}(l)];
                    plot_z_nH = [plot_z_nH, subvolume_xyz{1,3}(m)];
                    plot_nH = [plot_nH, subvolume_nH(k,l,m)];
                end
                if subvolume_incluster(k, l, m) > 0
                    plot_x_incluster = [plot_x_incluster, subvolume_xyz{1,1}(k)];
                    plot_y_incluster = [plot_y_incluster, subvolume_xyz{1,2}(l)];
                    plot_z_incluster = [plot_z_incluster, subvolume_xyz{1,3}(m)];
                    plot_incluster = [plot_incluster, subvolume_incluster(k,l,m)];
                end
                if removedH(k, l, m) > 0
                    plot_x_removedH = [plot_x_removedH, subvolume_xyz{1,1}(k)];
                    plot_y_removedH = [plot_y_removedH, subvolume_xyz{1,2}(l)];
                    plot_z_removedH = [plot_z_removedH, subvolume_xyz{1,3}(m)];
                    plot_removedH = [plot_removedH, removedH(k,l,m)];
                end
            end
        end
    end
    
    if plotflag == 1
             
        figure(2)
        scatter3(plot_x_nH, plot_y_nH, plot_z_nH, plot_nH)
        set(gcf, 'Visible', 'off')
        axis([0 N 0 N 0 M])
        %shg
        
        filename1 = sprintf('C:\\Users\\lk3618\\OneDrive - Imperial College London\\Sims\\ISFETSim\\SimulatingLAMPDiffusionGibsonBruck\\Results\\GibsonBruck_30_04_2020\\ProtonLocations_t%d.fig', int16(t(i)));
        savefig(filename1)
       
        figure(4)
        scatter3(plot_x_incluster, plot_y_incluster, plot_z_incluster, plot_incluster)
        set(gcf, 'Visible', 'off')
        axis([0 N 0 N 0 M])
        %shg
        
        filename3 = sprintf('C:\\Users\\lk3618\\OneDrive - Imperial College London\\Sims\\ISFETSim\\SimulatingLAMPDiffusionGibsonBruck\\Results\\GibsonBruck_30_04_2020\\Clusters_t%d.fig', int16(t(i)));
        savefig(filename3)
        
        figure(5)
        scatter3(plot_x_removedH, plot_y_removedH, plot_z_removedH, plot_removedH)
        set(gcf, 'Visible', 'off')
        axis([0 N 0 N 0 M])
        %shg
        
        
        sensornH_cell = mat2cell(subvolume_nH(:,:,1), xstep*ones(1, N), ystep*ones(1, N));
        sensornH = zeros(N-1, N-1);
        for k = 1:size(sensornH_cell, 1)
            for l = 1:size(sensornH_cell, 2)
                sensornH(k, l) = sum(sum(sensornH_cell{k, l}));
            end
        end
        
        figure(3)
        surface(sensor_meshy, sensor_meshx, sensornH)
        set(gcf, 'Visible', 'off')
        [value_max, ~] = max(sensornH(:));
        if value_max == 0
            value_max = 1
        end
        axis([0 N 0 N 0 value_max])
        %shg      
        
        filename2 = sprintf('C:\\Users\\lk3618\\OneDrive - Imperial College London\\Sims\\ISFETSim\\SimulatingLAMPDiffusionGibsonBruck\\Results\\GibsonBruck_30_04_2020\\Signal_t%d.fig', int16(t(i)));
        savefig(filename2)
    end
    
    %{
    [X,Y] = meshgrid(1:N);
    h = surf(X,Y,ones(N), 'FaceColor',[204, 255, 230]/255);
    %set(h,'facealpha',0.8);
    box on;
    
    drawnow limitrate;
    
    % Update sensor plot
    Z = sumSensor2(clusters, N);
    
    % Update count of real total # of protons
    tot_protons = sum(Z(:)) + persisted_protons;
    
    % Keep track of sensor array along time
    track(i,:,:) = Z;
    
    100*i/length(t);
    
    mystr = {'Release of Protons in 3D Space',['t = ' num2str(t(i)/60) ' minutes']};
    title(mystr, 'FontSize', 20);
    
    %mymov(i) = getframe(f);
    %}
end

toc;
% 
% %%
% v = VideoWriter('multicluster.avi');
% open(v);
% writeVideo(v, mymov);
% close(v);