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

%Record I_post
for n=1:length(event)-1
    event(n).I_post=event(n+1).startTime-event(n).endTime;
end

%Note if events are periodic. An event is periodic if Ipre and Ipost agree
%to within 20 percent 
for n=1:length(event)
    if abs(event(n).I_pre-event(n).I_post)/event(n).I_pre < 0.2
        event(n).periodic=true;
    end
end
%Count periodic events
numPeriodicEvents=0;
for n=1:length(event)
    if event(n).periodic
        numPeriodicEvents=numPeriodicEvents+1;
    end
end


%{
%Make a plot of packing angle vs. time
packingAngles=zeros(N,1);
for n=1:N
    packingAngles(n)=packing(n).packingAngleDegrees;
end
figure(1)
plot(1:N,packingAngles)
title("Packing Angle vs. Time")
xlabel("Time(s)")
ylabel("Packing Angle(Degrees)")
%}


%Make a histogram of maxima and mima on the same axes
maximaAngles=zeros(1,length(event));
minimaAngles=zeros(1,length(event));
for n=1:length(event)
    maximaAngles(n)=event(n).maxAngle;
    minimaAngles(n)=event(n).minAngle;
end
meanMinAngle=mean(minimaAngles);
sdMinAngle=std(minimaAngles);
skewnessMinAngles=skewness(minimaAngles);
kurtosisMinAngles=kurtosis(minimaAngles);
meanMaxAngle=mean(maximaAngles);
sdMaxAngle=std(maximaAngles);
skewnessMaxAngles=skewness(maximaAngles);
kurtosisMaxAngles=kurtosis(maximaAngles);
figure(2)
minHist=histogram(minimaAngles,'facealpha',.7,'edgecolor','none');
minHist;
hold on
maxHist=histogram(maximaAngles,'facealpha',.7,'edgecolor','none');
maxHist;
title("Histogram of Packing Angles")
xlabel("Packing Angle")
ylabel("Number of Frames")
legend('Minima','Maxima','location','northwest')
minPeak=minHist.BinEdges(find(minHist.Values==max(minHist.Values))); %Identify the peak of each distribution
maxPeak=maxHist.BinEdges(find(maxHist.Values==max(maxHist.Values)));
str1=sprintf("Minima:\npeak=%d degrees\nmean=%.1f degrees\nsd=%.1f degrees\nskewness=%.1f\nkurtosis=%.1f", minPeak,meanMinAngle,sdMinAngle,skewnessMinAngles,kurtosisMinAngles);
str2=sprintf("Maxima:\npeak=%d degrees\nmean=%.1f degrees\nsd=%.1f degrees\nskewness=%.1f\nkurtosis=%.1f", maxPeak,meanMaxAngle,sdMaxAngle,skewnessMaxAngles,kurtosisMaxAngles);
text(22,100,str1)
text(40,100,str2)
hold off


%Cumulative histogram of minimum angles
figure(3)
histogram(minimaAngles,'Normalization','cdf')
title("Cumulative Distribution of Minimum Angles")
xlabel("Minimum Angles")
ylabel("Cumulative Density")


%Cumulative histogram of minimum angles
figure(4)
histogram(maximaAngles,'Normalization','cdf')
title("Cumulative Distribution of Maximum Angles")
xlabel("Maximum Angles")
ylabel("Cumulative Density")


%Make historgram of event size (measured in degrees)
eventSizes=zeros(length(event),1);
for n=1:length(event)
    eventSizes(n)=event(n).size;
end
meanEventSize=mean(eventSizes);
sdEventSize=std(eventSizes);
skewnessEventSize=skewness(eventSizes);
kurtosisEventSize=kurtosis(eventSizes);
figure(5)
eventSizeHist=histogram(eventSizes);
eventSizeHist;
title("Histogram of avalanche Sizes")
xlabel("Avalanche Size")
ylabel("Number of Avalanches")
eventSizePeak=eventSizeHist.BinEdges(find(eventSizeHist.Values==max(eventSizeHist.Values)));
str=sprintf("Minima:\npeak=%d degrees\nmean=%.1f degrees\nsd=%.1f degrees\nskewness=%.1f\nkurtosis=%.1f", eventSizePeak,meanEventSize,sdEventSize,skewnessEventSize,kurtosisEventSize);
text(14,60,str);


%Cumulative histogram of event size
figure(6)
histogram(eventSizes,'Normalization','cdf')
title("Cumulative Distribution of Avalanche Sizes")
xlabel("Avalanche Size")
ylabel("Cumulative Density")

%Semilog histogram of event size
figure(7)
histogram(eventSizes);
set(gca,'YScale','log');
title("Semi-log Histogram of Avalanche Size")
xlabel("Avalanche Size")
ylabel("Number of Avalanches")

%log-log histogram of event size
figure(8)
histogram(log10(eventSizes));
set(gca,'YScale','log')
title("log-log Histogram of Avalanche Size")
xlabel("log of Avalanche Size")
ylabel("Number of Avalanches")

%Make historgram of event sizes less than 4 degrees
smallEventSizes=[];
for n=1:length(eventSizes)
    if eventSizes(n)<4
    smallEventSizes=[smallEventSizes, eventSizes(n)];
    end
end
figure(9)
histogram(smallEventSizes)
title("Histogram of avalanche sizes less than 4 degrees")
xlabel("Avalanche Size")
ylabel("Number of Avalanches")


%Plot avalanche size vs. maxima angle
figure(10)
scatter(maximaAngles,eventSizes)
title("Avalanche Size vs. Maximum Angle")
xlabel("Packing Angle Before Avalanche")
ylabel("Avalanche Size")


%Store some information about each event (avalanche). An event is the time between
%maximum and successive minimum
function event = extractEventData(packing)
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
end