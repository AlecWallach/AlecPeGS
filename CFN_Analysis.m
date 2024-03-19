directory = 'D:\DCIM\211MSDCF\';
solved_directory = 'D:\DCIM\211MSDCF\SolvedNetworks';
files = dir(fullframe(solved_directory, 'DSC*postProcessingWorkspace.mat'));
nFrames = length(files);

packingStruct=load(fullfile(directory,'packingStruct.mat')).packing;


for frame = 1:nFrames

    % note the name of the file so we can search packingStruct for the
    % packing angle
    fileName = files(frame).name;
    
    % load post processing struct
    postProcessingFile = fullfile(solved_directory,files);
    postProcessingStruct = load(postProcessingFile);

    % load force adjacency matrix and binary adjacency matrix
    FAM = postProcessingStruct.W;
    BAM = postProcessingStruct.B;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analyze distribution of forces in contact network
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot mean coordination number as a function of critical angle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end