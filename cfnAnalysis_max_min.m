directory = "C:\Users\Squishfolk\Desktop\Alec\211MSDCF";
solved_directory = "C:\Users\Squishfolk\Desktop\Alec\211MSDCF\Solved Networks";
ImageFiles= dir(fullfile(directory,"DSC*.JPG"));
PostProcessingFiles = dir(fullfile(solved_directory, 'DSC*postProcessingWorkspace.mat'));
nFrames = length(PostProcessingFiles);

% load lists that tell us what the max/min frames are (fix packingStruct
% labeling to make this easier)
minFrames = load(fullfile(directory,'minFrames.mat'));
maxFrames = load(fullfile(directory,'maxFrames.mat'));

minPostProcessingFiles = [];
for n=1:length(PostProcessingFiles)
    for i=1:length(minFrames(1).minFrames)
        if [PostProcessingFiles(n).name(1:8),'.JPG'] == minFrames(1).minFrames(i)
            filePath = dir(fullfile(solved_directory,PostProcessingFiles(n).name));
            minPostProcessingFiles=[minPostProcessingFiles; filePath];
        end
    end
end
nMinFrames=length(minPostProcessingFiles);

maxPostProcessingFiles = [];
for n=1:length(PostProcessingFiles)
    for i=1:length(maxFrames(1).maxFrames)
        if [PostProcessingFiles(n).name(1:8),'.JPG'] == maxFrames(1).maxFrames(i)
            filePath = dir(fullfile(solved_directory,PostProcessingFiles(n).name));
            maxPostProcessingFiles=[maxPostProcessingFiles; filePath];
        end
    end
end
nMaxFrames=length(maxPostProcessingFiles);

packingStruct=load(fullfile(directory,'packingStruct.mat')).packing;

critical_angles=zeros(nMaxFrames,1);
avg_coordination_number_max=zeros(nMaxFrames,1);
avg_beta_max=zeros(nMaxFrames,1);
avg_alpha_max=zeros(nMaxFrames,1);
num_connected_components_max=zeros(nMaxFrames,1);
longest_connected_component_max=zeros(nMaxFrames,1);
three_cycle_density_max=zeros(nMaxFrames,1);
mean_three_cycle_stability_max=zeros(nMaxFrames,1);
avg_normal_force_max=zeros(nMaxFrames,1);
avg_tangential_force_max=zeros(nMaxFrames,1);

parfor frame = 1:nMaxFrames

    % note the name of the file so we can search packingStruct for the
    % packing angle
    ImageFileName = [maxPostProcessingFiles(frame).name(1:8),'.JPG'];
    
    % load post processing struct
    PostProcessingFileName=maxPostProcessingFiles(frame).name;
    postProcessingFile = fullfile(solved_directory,PostProcessingFileName);
    postProcessingStruct = load(postProcessingFile);

    % load force adjacency matrix and binary adjacency matrix
    FAM = postProcessingStruct.W;
    BAM = postProcessingStruct.B;

    %Determine the critical angle of the packing by searching the
    %packingStruct
    for i=1:length(packingStruct)
        if packingStruct(i).fileName(end-11:end) == ImageFileName
            critical_angles(frame) = packingStruct(i).packingAngle;
            break
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze distribution of forces in contact network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    %We only care about one triangle of the matrix since its symmetrical
    lowerFAM = tril(FAM);

    %Make a vector containing all of the non-zero forces
    forceVec=[];
    for i=1:size(lowerFAM,1)
        for j=1:size(lowerFAM,2)
            if lowerFAM(i,j)>0
                forceVec=[forceVec, lowerFAM(i,j)];
            end
        end
    end


    %Find the mean, standard deviation, kurtosis of force distribution
    meanForce = mean(forceVec);
    forceSD = std(forceVec);
    forceKurtosis = kurtosis(forceVec);

    %Make a histogram of forces and add info to the plot
    figure(1)
    histogram(forceVec,50)
    xlabel("Force (N)","FontSize",14)
    ylabel("Number of contacts","FontSize",14)
    forceStr=sprintf("mean force = %0.3f N",meanForce);
    sdStr=sprintf("standard deviation = %0.3f N",forceSD);
    kurtosisStr=sprintf("kurtosis = %0.3f",forceKurtosis);
    text(0.3,120,forceStr,"FontSize",14);
    text(0.3,115,sdStr,"FontSize",14);
    text(0.3,110,kurtosisStr,"FontSize",14);
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record average coordination number of the packing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    total_contacts=sum(sum(BAM));
    num_particles=length(postProcessingStruct.particle);
    avg_coordination_number_max(frame)=total_contacts/num_particles;
    
%{
    %Do this with G2 rather than forces
    total_contacts=0;
    particle=postProcessingStruct.particle;
    for i=1:length(particle)
        total_contacts = total_contacts + length(particle(i).neighbours)
    end
    num_particles=length(postProcessingStruct.particle);
    avg_coordination_number_max(frame)=total_contacts/num_particles;
%}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record average force orientation of the packing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    particle=postProcessingStruct.particle;
    N=length(particle);
    betas=[];
    alphas=[];
    for n=1:N
        betas=[betas, particle(n).betas];
        for i=1:length(particle(n).alphas) % for some reason some of the alphas are row vectors and some are column vectors
            alphas=[alphas,particle(n).alphas(i)];
        end
    end
    avg_beta_max(frame)=mean(betas);
    avg_alpha_max(frame)=mean(alphas);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record average normal and tangential force of the packing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    normal_forces=[];
    tangential_forces=[];
    for n=1:N
        normal_forces=[normal_forces, particle(n).forces.*cos(particle(n).alphas)];
        tangential_forces=[tangential_forces, particle(n).forces.*cos(particle(n).alphas)];
    end
    avg_normal_force(frame)=mean(normal_forces);
    avg_tangential_force(frame)=mean(tangential_forces);
    %}


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find number of connected components in the packing & longest
    % connected component
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bins = conncomp(graph(BAM));
    num_connected_components_max(frame)=max(bins);
    component_sizes=accumarray(bins',1);
    longest_connected_component_max(frame)=max(component_sizes);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find 3-cycle density and mean stability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = size(FAM,1);
    numThreeCycles = 0;
    cycle_stabilities=[];

    % Generate all possible combinations of three vertices
    combinations = nchoosek(1:j, 3);
    num_combinations = size(combinations, 1);

    %Check each combination for 3-cycle
    for i = 1:num_combinations
        vertices = combinations(i, :);
        % record the force on each edge
        force1 = FAM(vertices(1), vertices(2));
        force2 = FAM(vertices(2), vertices(3));
        force3 = FAM(vertices(3), vertices(1));
        if force1~=0 && force2~=0 && force3~=0 
            % Found a 3-cycle
            numThreeCycles = numThreeCycles + 1;
            avg_force=(force1+force2+force3)/3;
            cycle_stability= (force1*force2*force3)/avg_force;
            cycle_stabilities=[cycle_stabilities,cycle_stability];
        end
    end

    three_cycle_density_max(frame)=numThreeCycles/N;
    mean_three_cycle_stability_max(frame)=mean(cycle_stabilities);


end



repose_angles=zeros(nMinFrames,1);
avg_coordination_number_min=zeros(nMaxFrames,1);
avg_beta_min=zeros(nMinFrames,1);
avg_alpha_min=zeros(nMinFrames);
num_connected_components_min=zeros(nMinFrames,1);
longest_connected_component_min=zeros(nMinFrames,1);
three_cycle_density_min=zeros(nMinFrames,1);
mean_three_cycle_stability_min=zeros(nMinFrames,1);
avg_normal_force_min=zeros(nMinFrames,1);
avg_tangential_force_min=zeros(nMinFrames,1);

parfor frame = 1:nMinFrames

    % note the name of the file so we can search packingStruct for the
    % packing angle
    ImageFileName = [minPostProcessingFiles(frame).name(1:8),'.JPG'];
    
    % load post processing struct
    PostProcessingFileName=minPostProcessingFiles(frame).name;
    postProcessingFile = fullfile(solved_directory,PostProcessingFileName);
    postProcessingStruct = load(postProcessingFile);

    % load force adjacency matrix and binary adjacency matrix
    FAM = postProcessingStruct.W;
    BAM = postProcessingStruct.B;

    %Determine the critical angle of the packing by searching the
    %packingStruct
    for i=1:length(packingStruct)
        if packingStruct(i).fileName(end-11:end) == ImageFileName
            repose_angles(frame) = packingStruct(i).packingAngle;
            break
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze distribution of forces in contact network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    %We only care about one triangle of the matrix since its symmetrical
    lowerFAM = tril(FAM);

    %Make a vector containing all of the non-zero forces
    forceVec=[];
    for i=1:size(lowerFAM,1)
        for j=1:size(lowerFAM,2)
            if lowerFAM(i,j)>0
                forceVec=[forceVec, lowerFAM(i,j)];
            end
        end
    end


    %Find the mean, standard deviation, kurtosis of force distribution
    meanForce = mean(forceVec);
    forceSD = std(forceVec);
    forceKurtosis = kurtosis(forceVec);

    %Make a histogram of forces and add info to the plot
    figure(1)
    histogram(forceVec,50)
    xlabel("Force (N)","FontSize",14)
    ylabel("Number of contacts","FontSize",14)
    forceStr=sprintf("mean force = %0.3f N",meanForce);
    sdStr=sprintf("standard deviation = %0.3f N",forceSD);
    kurtosisStr=sprintf("kurtosis = %0.3f",forceKurtosis);
    text(0.3,120,forceStr,"FontSize",14);
    text(0.3,115,sdStr,"FontSize",14);
    text(0.3,110,kurtosisStr,"FontSize",14);
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record average coordination number of the packing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    total_contacts=sum(sum(BAM));
    num_particles=length(postProcessingStruct.particle);
    avg_coordination_number_min(frame)=total_contacts/num_particles;
    
%{
    %Do this with G2 rather than forces
    total_contacts=0;
    particle=postProcessingStruct.particle;
    for i=1:length(particle)
        total_contacts = total_contacts + length(particle(i).neighbours)
    end
    num_particles=length(postProcessingStruct.particle);
    avg_coordination_number_min(frame)=total_contacts/num_particles;
%}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record average force orientation of the packing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    particle=postProcessingStruct.particle;
    N=length(particle);
    betas=[];
    alphas=[];
    for n=1:N
        betas=[betas, particle(n).betas];
        for i=1:length(particle(n).alphas) % for some reason some of the alphas are row vectors and some are column vectors
            alphas=[alphas,particle(n).alphas(i)];
        end
    end
    avg_beta_min(frame)=mean(betas);
    avg_alpha_min(frame)=mean(alphas);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record average normal and tangential force of the packing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    normal_forces=[];
    tangential_forces=[];
    for n=1:N
        normal_forces=[normal_forces, particle(n).forces.*cos(particle(n).alphas)];
        tangential_forces=[tangential_forces, particle(n).forces.*cos(particle(n).alphas)];
    end
    avg_normal_force(frame)=mean(normal_forces);
    avg_tangential_force(frame)=mean(tangential_forces);
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find number of connected components in the packing & longest
    % connected component
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bins = conncomp(graph(BAM));
    num_connected_components_min(frame)=max(bins);
    component_sizes=accumarray(bins',1);
    longest_connected_component_min(frame)=max(component_sizes);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find 3-cycle density and mean stability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = size(FAM,1);
    numThreeCycles = 0;
    cycle_stabilities=[];

    % Generate all possible combinations of three vertices
    combinations = nchoosek(1:j, 3);
    num_combinations = size(combinations, 1);

    %Check each combination for 3-cycle
    for i = 1:num_combinations
        vertices = combinations(i, :);
        % record the force on each edge
        force1 = FAM(vertices(1), vertices(2));
        force2 = FAM(vertices(2), vertices(3));
        force3 = FAM(vertices(3), vertices(1));
        if force1~=0 && force2~=0 && force3~=0 
            % Found a 3-cycle
            numThreeCycles = numThreeCycles + 1;
            avg_force=(force1+force2+force3)/3;
            cycle_stability= (force1*force2*force3)/avg_force;
            cycle_stabilities=[cycle_stabilities,cycle_stability];
        end
    end

    three_cycle_density_min(frame)=numThreeCycles/N;
    mean_three_cycle_stability_min(frame)=mean(cycle_stabilities);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordination number plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot mean coordination number vs. critical angle
figure(1)
subplot(1,2,1);
scatter(critical_angles,avg_coordination_number_max);
hold on
scatter(repose_angles,avg_coordination_number_min);
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
xline(mean(repose_angles),'--','\color{red}AR', 'FontSize', 14)
xlabel('$\theta$','Interpreter','latex', 'FontSize', 18)
ylabel('Mean Coordination Number', 'FontSize', 14)
box on

% Plot the centroid for the maxima and minima
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles) / sqrt(length(critical_angles));
SEM_y = std(avg_coordination_number_max) / sqrt(length(avg_coordination_number_max));
scatter(mean(critical_angles),mean(avg_coordination_number_max),'filled', 'blue', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles), mean(avg_coordination_number_max), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

SEM_x = std(repose_angles) / sqrt(length(repose_angles));
SEM_y = std(avg_coordination_number_min) / sqrt(length(avg_coordination_number_min));
scatter(mean(repose_angles),mean(avg_coordination_number_min),'filled', 'red', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(repose_angles), mean(avg_coordination_number_min), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

hold off

% Create the difference plot
subplot(1,2,2)
coordination_number_difference=[];
for i=2:length(avg_coordination_number_max)
    max_coordination_number=avg_coordination_number_max(i);
    min_coordination_number=avg_coordination_number_min(i-1);
    coordination_number_difference=[coordination_number_difference, max_coordination_number-min_coordination_number];
end
scatter(critical_angles(2:end),coordination_number_difference)
xlabel('MSA', 'FontSize', 14)
ylabel('Change in Mean Coordination Number', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
box on

% Plot the centroid of the difference plot
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles(2:end)) / sqrt(length(critical_angles(2:end)));
SEM_y = std(coordination_number_difference) / sqrt(length(coordination_number_difference));
hold on
scatter(mean(critical_angles(2:end)),mean(coordination_number_difference),'filled', 'black', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles(2:end)), mean(coordination_number_difference), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off

sgtitle('Coordination Number','FontSize', 18)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beta angle plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot mean beta angle vs. critical angle
figure(2)
subplot(1,2,1);
scatter(critical_angles,avg_beta_max);
hold on
scatter(repose_angles,avg_beta_min);
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
xline(mean(repose_angles),'--','\color{red}AR', 'FontSize', 14)
xlabel('$\theta$','Interpreter','latex', 'FontSize', 18)
ylabel('Mean Beta Angle', 'FontSize', 14)
box on

% Plot the centroid for the maxima and minima
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles) / sqrt(length(critical_angles));
SEM_y = std(avg_beta_max) / sqrt(length(avg_beta_max));
scatter(mean(critical_angles),mean(avg_beta_max),'filled', 'blue', 'SizeData', 75)
% Add error bars to the centroid
errorbar(mean(critical_angles), mean(avg_beta_max), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

SEM_x = std(repose_angles) / sqrt(length(repose_angles));
SEM_y = std(avg_beta_min) / sqrt(length(avg_beta_min));
scatter(mean(repose_angles),mean(avg_beta_min),'filled', 'red', 'SizeData', 75)
% Add error bars to the centroid
errorbar(mean(repose_angles), mean(avg_beta_min), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

hold off

% Create the difference plot
subplot(1,2,2)
beta_angle_difference=[];
for i=2:length(avg_beta_max)
    max_beta_angle=avg_beta_max(i);
    min_beta_angle=avg_beta_min(i-1);
    beta_angle_difference=[beta_angle_difference, max_beta_angle-min_beta_angle];
end
scatter(critical_angles(2:end),beta_angle_difference)
xlabel('MSA', 'FontSize', 14)
ylabel('Change in Mean Beta Angle', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
box on

% Plot the centroid of the difference plot
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles(2:end)) / sqrt(length(critical_angles(2:end)));
SEM_y = std(beta_angle_difference) / sqrt(length(beta_angle_difference));
hold on
scatter(mean(critical_angles(2:end)),mean(beta_angle_difference),'filled', 'black', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles(2:end)), mean(beta_angle_difference), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off

sgtitle('Beta Angle','FontSize', 18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connected Components plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot number of connected components vs.angle
figure(3)
subplot(2,2,1);
scatter(critical_angles,num_connected_components_max);
hold on
scatter(repose_angles,num_connected_components_min);
ylabel('Number of Connected Components', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
xline(mean(repose_angles),'--','\color{red}AR', 'FontSize', 14)
box on

% Plot the centroid for the maxima and minima
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles) / sqrt(length(critical_angles));
SEM_y = std(num_connected_components_max) / sqrt(length(num_connected_components_max));
scatter(mean(critical_angles),mean(num_connected_components_max),'filled', 'blue', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles), mean(num_connected_components_max), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

disp(max(critical_angles))
disp(SEM_x)

SEM_x = std(repose_angles) / sqrt(length(repose_angles));
SEM_y = std(num_connected_components_min) / sqrt(length(num_connected_components_min));
scatter(mean(repose_angles),mean(num_connected_components_min),'filled', 'red', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(repose_angles), mean(num_connected_components_max), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off

% Plot longest connected component vs.angle
subplot(2,2,3);
scatter(critical_angles,longest_connected_component_max);
hold on
scatter(repose_angles,longest_connected_component_min);
xlabel('$\theta$','Interpreter','latex', 'FontSize', 18)
ylabel('Length of Longest Connected Component', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
xline(mean(repose_angles),'--','\color{red}AR', 'FontSize', 14)

box on

% Plot the centroid for the maxima and minima
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles) / sqrt(length(critical_angles));
SEM_y = std(longest_connected_component_max) / sqrt(length(longest_connected_component_max));
scatter(mean(critical_angles),mean(longest_connected_component_max),'filled', 'blue', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles), mean(longest_connected_component_max), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

SEM_x = std(repose_angles) / sqrt(length(repose_angles));
SEM_y = std(longest_connected_component_min) / sqrt(length(longest_connected_component_min));
scatter(mean(repose_angles),mean(longest_connected_component_min),'filled', 'red', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(repose_angles), mean(longest_connected_component_min), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off


% Plot difference between number of connected components for frame before avalanche n
% and coordination number for frame after avalanche n-1 vs. critical angle
subplot(2,2,2)
num_connected_components_difference=[];
for i=2:length(num_connected_components_max)
    max_num_connected_components=num_connected_components_max(i);
    min_num_connected_components=num_connected_components_min(i-1);
    num_connected_components_difference=[num_connected_components_difference, max_num_connected_components-min_num_connected_components];
end
scatter(critical_angles(2:end),num_connected_components_difference)
ylabel('Change in Number of Connected Components', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)

box on

% Plot the centroid of the difference plot
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles(2:end)) / sqrt(length(critical_angles(2:end)));
SEM_y = std(num_connected_components_difference) / sqrt(length(num_connected_components_difference));
hold on
scatter(mean(critical_angles(2:end)),mean(num_connected_components_difference),'filled', 'black', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles(2:end)), mean(num_connected_components_difference), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off

% Plot difference between coordination number for frame before avalanche n
% and coordination number for frame after avalanche n-1 vs. critical angle
subplot(2,2,4)
longest_connected_component_difference=[];
for i=2:length(longest_connected_component_max)
    max_longest_connected_component=longest_connected_component_max(i);
    min_longest_connected_component=longest_connected_component_min(i-1);
    longest_connected_component_difference=[longest_connected_component_difference, max_longest_connected_component-min_longest_connected_component];
end
scatter(critical_angles(2:end),longest_connected_component_difference)
xlabel('MSA', 'FontSize', 14)
ylabel('Change in Length of Longest Connected Component', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)

box on

% Plot the centroid of the difference plot
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles(2:end)) / sqrt(length(critical_angles(2:end)));
SEM_y = std(longest_connected_component_difference) / sqrt(length(longest_connected_component_difference));
hold on
scatter(mean(critical_angles(2:end)),mean(longest_connected_component_difference),'filled', 'black', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles(2:end)), mean(longest_connected_component_difference), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off


sgtitle('Connected Components','FontSize',18)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-cycle plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot 3-cycle density vs.angle
figure(4)
subplot(2,2,1);
scatter(critical_angles,three_cycle_density_max);
hold on
scatter(repose_angles,three_cycle_density_min);
ylabel('3-Cycle Density', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
xline(mean(repose_angles),'--','\color{red}AR', 'FontSize', 14)

box on

% Plot the centroid for the maxima and minima
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles) / sqrt(length(critical_angles));
SEM_y = std(three_cycle_density_max) / sqrt(length(three_cycle_density_max));
scatter(mean(critical_angles),mean(three_cycle_density_max),'filled', 'blue', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles), mean(three_cycle_density_max), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

disp(max(critical_angles))
disp(SEM_x)

SEM_x = std(repose_angles) / sqrt(length(repose_angles));
SEM_y = std(three_cycle_density_min) / sqrt(length(three_cycle_density_min));
scatter(mean(repose_angles),mean(three_cycle_density_min),'filled', 'red', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(repose_angles), mean(three_cycle_density_min), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off

% Plot longest connected component vs.angle
subplot(2,2,3);
scatter(critical_angles,mean_three_cycle_stability_max);
hold on
scatter(repose_angles,mean_three_cycle_stability_min);
xlabel('$\theta$','Interpreter','latex', 'FontSize', 18)
ylabel('Mean 3-cycle Stability', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)
xline(mean(repose_angles),'--','\color{red}AR', 'FontSize', 14)

box on

% Plot the centroid for the maxima and minima
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles) / sqrt(length(critical_angles));
SEM_y = std(mean_three_cycle_stability_max) / sqrt(length(mean_three_cycle_stability_max));
scatter(mean(critical_angles),mean(mean_three_cycle_stability_max),'filled', 'blue', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles), mean(mean_three_cycle_stability_max), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);

SEM_x = std(repose_angles) / sqrt(length(repose_angles));
SEM_y = std(mean_three_cycle_stability_min) / sqrt(length(mean_three_cycle_stability_min));
scatter(mean(repose_angles),mean(mean_three_cycle_stability_min),'filled', 'red', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(repose_angles), mean(mean_three_cycle_stability_min), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off


% Plot difference between number of connected components for frame before avalanche n
% and coordination number for frame after avalanche n-1 vs. critical angle
subplot(2,2,2)
three_cycle_density_difference=[];
for i=2:length(three_cycle_density_max)
    max_three_cycle_density=three_cycle_density_max(i);
    min_three_cycle_density=three_cycle_density_min(i-1);
    three_cycle_density_difference=[three_cycle_density_difference, max_three_cycle_density-min_three_cycle_density];
end
scatter(critical_angles(2:end),three_cycle_density_difference)
ylabel('Change in 3-Cycle Density', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)

box on

% Plot the centroid of the difference plot
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles(2:end)) / sqrt(length(critical_angles(2:end)));
SEM_y = std(three_cycle_density_difference) / sqrt(length(three_cycle_density_difference));
hold on
scatter(mean(critical_angles(2:end)),mean(three_cycle_density_difference),'filled', 'black', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles(2:end)), mean(three_cycle_density_difference), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off

% Plot difference between 3-cycle stability for frame before avalanche n
% and for frame after avalanche n-1 vs. critical angle
subplot(2,2,4)
three_cycle_stability_difference=[];
for i=2:length(mean_three_cycle_stability_max)
    max_three_cycle_stability=mean_three_cycle_stability_max(i);
    min_three_cycle_stability=mean_three_cycle_stability_min(i-1);
    three_cycle_stability_difference=[three_cycle_stability_difference, max_three_cycle_stability-min_three_cycle_stability];
end
scatter(critical_angles(2:end),three_cycle_stability_difference)
xlabel('MSA', 'FontSize', 14)
ylabel('Change in Mean 3-cycle Stability', 'FontSize', 14)
xline(mean(critical_angles),'--','\color{blue}MSA', 'FontSize', 14)

box on

% Plot the centroid of the difference plot
% Calculate standard error of the mean for x and y coordinates of the centroid
SEM_x = std(critical_angles(2:end)) / sqrt(length(critical_angles(2:end)));
SEM_y = std(three_cycle_stability_difference) / sqrt(length(three_cycle_stability_difference));
hold on
scatter(mean(critical_angles(2:end)),mean(three_cycle_stability_difference),'filled', 'black', 'SizeData', 60)
% Add error bars to the centroid
errorbar(mean(critical_angles(2:end)), mean(three_cycle_stability_difference), ...
    SEM_y, SEM_y, SEM_x, SEM_x, 'k', 'LineStyle', 'none', 'LineWidth', 1);
hold off


sgtitle('3-Cycles','FontSize',18)



%{
% Plot mean beta angle vs.critical angle
figure(2)
scatter(critical_angles,avg_beta_max);
xlabel('Critical Angle')
ylabel('Mean Beta Angle')
title('Mean Beta Angle vs. Critical Angle')

% Plot mean alpha angle vs.critical angle
figure(3)
scatter(critical_angles,avg_alpha_max);
xlabel('Critical Angle')
ylabel('Mean Alpha Angle')
title('Mean Alpha Angle vs. Critical Angle')

% Plot number of connected components vs.critical angle
figure(4)
scatter(critical_angles,num_connected_components_max);
xlabel('Critical Angle')
ylabel('Number of Connected Components')
title('Number of Connected Components vs. Critical Angle')

% Plot longest connected components vs.critical angle
figure(5)
scatter(critical_angles,longest_connected_component_max);
xlabel('Critical Angle')
ylabel('Longest Connected Components')
title('NLongest Connected Component vs. Critical Angle')


% Plot 3-cycle density vs.critical angle
figure(6)
scatter(critical_angles,three_cycle_density_max);
xlabel('Critical Angle')
ylabel('3-cycle density')
title('3-cycle density vs. Critical Angle')

% Plot mean 3-cycle stability vs.critical angle
figure(7)
scatter(critical_angles,mean_three_cycle_stability_max);
xlabel('Critical Angle')
ylabel('mean 3-cycle stability')
title('mean 3-cycle stability vs. Critical Angle')

% Plot mean coordination number vs. angle with maxs and mins in different
% colors
figure(8)
scatter(critical_angles,avg_coordination_number_max);
hold on
scatter(repose_angles,avg_coordination_number_min);
xlabel('Critical Angle')
ylabel('Mean Coordination Number')
title('Mean Coordination Number vs. Critical Angle')
hold off

% Plot mean beta vs. angle with maxs and mins in different
% colors
figure(9)
scatter(critical_angles,avg_beta_max);
hold on
scatter(repose_angles,avg_beta_min);
xlabel('Critical Angle')
ylabel('Mean Beta Angle')
title('Mean Beta Angle vs. Critical Angle')
hold off

% Plot mean alpha angle vs. angle with maxs and mins in different
% colors
figure(10)
scatter(critical_angles,avg_alpha_max);
hold on
scatter(repose_angles,avg_alpha_min);
xlabel('Critical Angle')
ylabel('Mean Alpha Angle')
title('Mean Alpha Angle vs. Critical Angle')
hold off

% Plot number of connected components vs. angle with maxs and mins in different
% colors
figure(11)
scatter(critical_angles,num_connected_components_max);
hold on
scatter(repose_angles,num_connected_components_min);
xlabel('Critical Angle')
ylabel('Number of Connected Components')
title('Number of Connected Components vs. Critical Angle')
hold off

% Plot longest connected components vs. angle with maxs and mins in different
% colors
figure(12)
scatter(critical_angles,longest_connected_component_max);
hold on
scatter(repose_angles,longest_connected_component_min);
xlabel('Critical Angle')
ylabel('Longest Connected Component')
title('Longest Connected Component vs. Critical Angle')
hold off

% Plot mean alpha angle vs. angle with maxs and mins in different
% colors
figure(13)
scatter(critical_angles,avg_alpha_max);
hold on
scatter(repose_angles,avg_alpha_min);
xlabel('Critical Angle')
ylabel('Mean Alpha Angle')
title('Mean Alpha Angle vs. Critical Angle')
hold off


% Plot 3-cycle density vs.critical angle with maxs and mins in different
% colors
figure(14)
scatter(critical_angles,three_cycle_density_max);
hold on
scatter(repose_angles,three_cycle_density_min);
xlabel('Critical Angle')
ylabel('3-cycle density')
title('3-cycle density vs. Critical Angle')
hold off

% Plot mean 3-cycle stability vs.critical angle with maxs and mins in different
% colors
figure(15)
scatter(critical_angles,mean_three_cycle_stability_max);
hold on
scatter(repose_angles,mean_three_cycle_stability_min);
xlabel('Critical Angle')
ylabel('mean 3-cycle stability')
title('mean 3-cycle stability vs. Critical Angle')
hold off


% Plot difference between coordination number for frame before avalanche n
% and coordination number for frame after avalanche n-1 vs. critical angle
coordination_number_difference=[];
for i=2:length(avg_coordination_number_max)
    max_coordination_number=avg_coordination_number_max(i);
    min_coordination_number=avg_coordination_number_min(i-1);
    coordination_number_difference=[coordination_number_difference, max_coordination_number-min_coordination_number];
end
figure(16)
scatter(critical_angles(2:end),coordination_number_difference)
xlabel('Critical Angle of Avalanche n')
ylabel('Difference in Mean Coordination Number')
title('Coordination Number Difference vs. Critical Angle')

% Plot difference between beta angle for frame before avalanche n
% and beta angle for frame after avalanche n-1 vs. critical angle
beta_angle_difference=[];
for i=2:length(avg_beta_max)
    max_beta_angle=avg_beta_max(i);
    min_beta_angle=avg_beta_min(i-1);
    beta_angle_difference=[beta_angle_difference, max_beta_angle-min_beta_angle];
end
figure(17)
scatter(critical_angles(2:end),beta_angle_difference)
xlabel('Critical Angle of Avalanche n')
ylabel('Difference in Mean Beta Angle')
title('Mean Beta Angle Difference vs. Critical Angle')

% Plot difference between alpha angle for frame before avalanche n
% and alpha for frame after avalanche n-1 vs. critical angle
alpha_angle_difference=[];
for i=2:length(avg_alpha_max)
    max_alpha_angle=avg_alpha_max(i);
    min_alpha_angle=avg_alpha_min(i-1);
    alpha_angle_difference=[alpha_angle_difference, max_alpha_angle-min_alpha_angle];
end
figure(18)
scatter(critical_angles(2:end),alpha_angle_difference)
xlabel('Critical Angle of Avalanche n')
ylabel('Difference in Mean alpha Angle')
title('Mean alpha Angle Difference vs. Critical Angle')

% Plot difference between number of connected components for frame before avalanche n
% and number of connected components for frame after avalanche n-1 vs. critical angle
num_connected_components_difference=[];
for i=2:length(num_connected_components_max)
    max_num_connected_components=num_connected_components_max(i);
    min_num_connected_components=num_connected_components_min(i-1);
    num_connected_components_difference=[num_connected_components_difference, max_num_connected_components-min_num_connected_components];
end
figure(19)
scatter(critical_angles(2:end),num_connected_components_difference)
xlabel('Critical Angle of Avalanche n')
ylabel('Difference in Number of Connected Components')
title('Number of Connected Components Difference vs. Critical Angle')

% Plot difference between length of longest connected component for frame before avalanche n
% and longest connected component for frame after avalanche n-1 vs. critical angle
longest_connected_component_difference=[];
for i=2:length(longest_connected_component_max)
    max_longest_connected_component=longest_connected_component_max(i);
    min_longest_connected_component=longest_connected_component_min(i-1);
    longest_connected_component_difference=[longest_connected_component_difference, max_longest_connected_component-min_longest_connected_component];
end
figure(20)
scatter(critical_angles(2:end),longest_connected_component_difference)
xlabel('Critical Angle of Avalanche n')
ylabel('Difference in Length of Longest Connected Component')
title('Length of Longest Connected Component Difference vs. Critical Angle')

% Plot difference between 3-cycle density for frame before avalanche n
% and for frame after avalanche n-1 vs. critical angle
three_cycle_density_difference=[];
for i=2:length(three_cycle_density_max)
    max_three_cycle_density=three_cycle_density_max(i);
    min_three_cycle_density=three_cycle_density_min(i-1);
    three_cycle_density_difference=[three_cycle_density_difference, max_three_cycle_density-min_three_cycle_density];
end
figure(21)
scatter(critical_angles(2:end),three_cycle_density_difference)
xlabel('Critical Angle of Avalanche n')
ylabel('Difference in 3-cycle density')
title('3-cycle density Difference vs. Critical Angle')

% Plot difference between 3-cycle stability for frame before avalanche n
% and for frame after avalanche n-1 vs. critical angle
three_cycle_stability_difference=[];
for i=2:length(mean_three_cycle_stability_max)
    max_three_cycle_stability=mean_three_cycle_stability_max(i);
    min_three_cycle_stability=mean_three_cycle_stability_min(i-1);
    three_cycle_stability_difference=[three_cycle_stability_difference, max_three_cycle_stability-min_three_cycle_stability];
end
figure(22)
scatter(critical_angles(2:end),three_cycle_stability_difference)
xlabel('Critical Angle of Avalanche n')
ylabel('Difference in mean 3-cycle stability')
title('Mean 3-cycle stability Difference vs. Critical Angle')
%}
