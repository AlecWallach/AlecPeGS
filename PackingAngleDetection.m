%This program takes raw and preprocessed packing images as inputs and saves a
%file containing the slope of each packing. The preprocessing should create
%as much contrast as possible between the particles and their background 
%(e.g. green channel, threshold, blur). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directory = '/Volumes/SD/DCIM/211MSDCF/';
preProDirectory='/Volumes/SD/DCIM/211MSDCF/PreProImages/';
files = dir(fullfile(directory, 'DSC*.JPG'));
preProFiles=dir(fullfile(preProDirectory,'DSC*.jpg'));
nFrames = length(files); %how many files are we processing ?

%Do we want to make angle vs. time animation
makeAngleTimeAnimation = false;

%How many frames were taken per second
fps = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End User Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear packing;

%Create a structure array to store information about each frame
packing(1:nFrames) = struct('fileName',[],'frameNumber',[],'time',[],'slope',[],'intercept',[],'packingAngle',[],'packingAngleDegrees',[],'maximum',[],'minimum',[]);

for frame = 1:nFrames %Loops for total number of images

    %This just outputs the progress of the program
    if mod(frame,50)==0
        disp(frame)
    end

    %Access the relevant image
    imageFile = [directory,files(frame).name]; %input filename
    imageFilePrePro = [preProDirectory,preProFiles(frame).name];
    img = imread(imageFile);
    preProImg = imread(imageFilePrePro);

    %Store some information about the packing
    packing(frame).fileName=imageFile;
    packing(frame).frameNumber=frame;
    packing(frame).time=frame*fps;
    
    Rimg = preProImg;
    
    %Crop Image to just the portion of the drum that contains particles
    %by drawing circle on first frame in dataset
    if frame==1
    imshow(Rimg)
    drumCrop = drawcircle;
    drumCenter=round(drumCrop.Center);
    drumRadius=drumCrop.Radius;
    end
    Rimg=imadjust(Rimg,[0.1, 1]);
    Rimg=imbinarize(Rimg);


    %Look at the middle two thirds of the x coordinates of the drum.
    %Divide this area into ten vertical strips of equal width.
    %Divide those vertical strips into 50 rectangles of equal height and
    %look at the intensity as a function of height. The point where the
    %intensity increases substantially is the top of the packing. If we fit
    %a line to all of these points, we should be able to approximate the
    %packing angle.
    xstart=round(drumCenter(1)-0.67*drumRadius);
    xend=round(drumCenter(1)+0.67*drumRadius);
    dx=round((xend-xstart)/10);
    Y=zeros(10,1);
    i=1;
    for x=xstart:dx:xend

        verticalStrip=Rimg(:,x:x+dx);
    
        dy=size(verticalStrip,1)/49;
        intensity_dA=zeros(50,1);
        n=1;
        for y = dy:dy:size(verticalStrip,1)-dy-1
            dA = verticalStrip(round(y):round(y+dy),:);
            intensity_dA(n) = sum(sum(dA));
            n=n+1;
        end

        for j=1:length(intensity_dA)-1
            if intensity_dA(j)>5000 && j>10
                yMax = (j)*dy;
                break
            end
        end
        Y(i)=yMax;
        i=i+1;
    end  

    %Fit a line to the points at the top of the packing
    X = xstart:dx:xend;
    P = polyfit(X,Y,1);

    %Record the slope, intercept, and packing angle
    packing(frame).slope=P(1);
    packing(frame).intercept=P(2);
    theta=atan(P(1));
    thetaDegrees=theta*180/pi;
    packing(frame).packingAngle=theta;
    packing(frame).packingAngleDegrees=theta*180/pi;



end

%Save packing structure array
save([directory,'packingStruct'],'packing');

