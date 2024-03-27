directory = "C:\Users\Squishfolk\Desktop\Alec\211MSDCF";
solved_directory = "C:\Users\Squishfolk\Desktop\Alec\211MSDCF\Solved Networks";
ImageFiles= dir(fullfile(directory,"DSC*.JPG"));
PostProcessingFiles = dir(fullfile(solved_directory, 'DSC*postProcessingWorkspace.mat'));
nFrames = length(files);

packingStruct=load(fullfile(directory,'packingStruct.mat')).packing;

critical_angles=zeros(nFrames,1);
avg_coordination_number=zeros(nFrames,1);

for frame = 1:nFrames

    % note the name of the file so we can search packingStruct for the
    % packing angle
    ImageFileName = ImageFiles(frame).name;
    
    % load post processing struct
    PostProcessingFileName=PostProcessingFiles(frame).name;
    postProcessingFile = fullfile(solved_directory,PostProcessingFileName);
    postProcessingStruct = load(postProcessingFile);

    % load force adjacency matrix and binary adjacency matrix
    FAM = postProcessingStruct.W;
    BAM = postProcessingStruct.B;

    %Determine the critical angle of the packing by searching the
    %packingStruct
    for i=1:length(packingStruct)
        if packingStruct(i).fileName(end-11:end) == ImageFileName
            critical_angles(frame) = packingStruct(i).packingAngleDegrees;
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
    avg_coordination_number(frame)=total_contacts/num_particles;

end

scatter(critical_angles,avg_coordination_number);