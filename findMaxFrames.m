directory = "C:\Users\Squishfolk\Desktop\Alec\211MSDCF";
fileName='eventStruct.mat';
file=load(fullfile(directory,fileName));
event=file.event;

maxFrames = strings(length(event),1);

for n=1:length(event)
    maxFrames(n,1) = event(n).startFrame(end-11:end);
end

disp(maxFrames)

%Save maxFrames
savePath = fullfile(directory, 'maxFrames');
% Save the file
save(savePath, 'maxFrames');