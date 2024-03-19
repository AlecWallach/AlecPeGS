%This program stores some information about each event (avalanche). An event is the time between
%maximum and successive minimum

%load the packing struct for the dataset of interest (the packing struct is 
%generated from PackingAngleDetection.m)
directory = '/Volumes/SD/DCIM/211MSDCF/';
fileName='packingStruct.mat';
file=load([directory,fileName]);
packing=file.packing;


%Store some information about each event (avalanche). An event is the time between
%maximum and successive minimum
N=length(packing);
eventNum=1;
n=1;
while n<N
    if packing(n+1).packingAngleDegrees < packing(n).packingAngleDegrees 
        
        %Find the frame where the slope stops decreasing i.e. the minumum
        j=n+1;
        if j<N-1
            while packing(j+1).packingAngleDegrees < packing(j).packingAngleDegrees
                j=j+1;
            end
        end

        packing(n).maximum = true;
        packing(j).minimum = true;
        
        event(eventNum).startFrame=packing(n).fileName;
        event(eventNum).startTime=packing(n).time;
        event(eventNum).endTime=packing(j).time;
        event(eventNum).maxAngle=packing(n).packingAngleDegrees;
        event(eventNum).minAngle=packing(j).packingAngleDegrees;
        event(eventNum).size=event(eventNum).maxAngle-event(eventNum).minAngle;
        if eventNum ~= 1
            event(eventNum).I_pre=event(eventNum).startTime-event(eventNum-1).endTime;
        end
        eventNum = eventNum+1;
        n=j;
    else
        n=n+1;
    end

end

%Save event structure array
save([directory,'eventStruct'],'event');