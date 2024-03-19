%Extract the maxima from the dataset

directory = '/Volumes/SD/DCIM/211MSDCF/';
packingFileName='packingStruct.mat';
eventFileName='eventStruct.mat';

clear packing;
clear event;

packingFile=load([directory,packingFileName]);
eventFile=load([directory,eventFileName]);
packing=packingFile.packing;
event=eventFile.event;

N=length(event);

maxFrames=[];

for n=1:N
    maxFrames=[maxFrames; event(n).startFrame];
end