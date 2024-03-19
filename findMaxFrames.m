directory = '/Volumes/SD/DCIM/211MSDCF/';
fileName='eventStruct.mat';
file=load([directory,fileName]);
event=file.event;

disp(event(1).startFrame)