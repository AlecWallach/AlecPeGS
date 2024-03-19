directory = '/Users/alecwallach/Desktop/Squishlab/Imaging/';
files = dir([directory, 'DSC_3852_preprocessing.mat']); 

nFrames = length(files); %how many files are we processing ?
for frame = 1:nFrames %loop over these frames 

    fileName = [directory,files(frame).name]; %which file/frame are we processing now ?
    load(fileName);

    

end