%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img = imread('/Users/alecwallach/Desktop/Squishlab/Imaging/DSC_3852.JPG');
img_prepro = imread('/Users/alecwallach/Desktop/Squishlab/Imaging/DSC_3852_Thresh_blur.jpg');


DS = 0.0025; % How much should we adjust sensitivity if wrong number of particles are found
RlargeH = [59 63]; %What radius (in pixels) range do we expect for the large discs?
RsmallH = [39 43]; %What radius (in pixels) range do we expect for the small discs?
SL = 0.975; %Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...
SS = 0.97; %Sensitivity of the Hough Transform disc detetcor, exact value is Voodo magic...

NsmallH = 262; %Number of small discs
NlargeH = 285; %Number of large discs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rimg = img_prepro(:,:,:);
Rimg = Rimg*0.4;

%Crop Image to just the portion of the drum that contains particles
imshow(Rimg)
dCrop = drawcircle; 
mask = createMask(dCrop);
Rimg = bsxfun(@times, Rimg, cast(mask, class(Rimg))); 
   
imshow(Rimg);
d = drawcircle;
RL = d.Radius;
RL = round(RL);

d = drawcircle;
RS = d.Radius;
RS = round(RS);

%Detect large and small particles with generous sensitivity
[centersL,radiiL,metricL] = imfindcircles(Rimg,[RL-2 RL+1],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.985,EdgeThreshold=0.025);
[centersS,radiiS,metricS] = imfindcircles(Rimg,[RS-2 RS+1],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.9775,EdgeThreshold=0.025);


%If two large particlces are overlapping by 20 percent of large radius,
%remove the one that is less circular (lower metric)
NL = length(centersL);
n=1;
while n<=NL
    j=n;
    while j<=NL
        if (sqrt(((centersL(n,1)-centersL(j,1)).^2)+((centersL(n,2)-centersL(j,2)).^2))<((2*RL)-(.20*RL))) && n~=j
            if(metricL(n)<metricL(j))
                centersL(n,:)=[];
                radiiL(n)=[];
                metricL(n)=[];
                NL=NL-1;
                n=n-1;
                break
            else
                centersL(j,:)=[];
                radiiL(j)=[];
                metricL(j)=[];
                NL=NL-1;
            end
         end
        j=j+1;
    end
    n=n+1;
end



%Remove weaker of two overlapping small particles
NS = length(centersS);
n=1;
while n<=NS
    j=n;
    while j<=NS
        if (sqrt(((centersS(n,1)-centersS(j,1))^2)+((centersS(n,2)-centersS(j,2))^2))<((2*RS)-(.20*RS))) && n~=j
            if(metricS(n)<metricS(j))
                centersS(n,:)=[];
                radiiS(n)=[];
                metricS(n)=[];
                NS=NS-1;
                n=n-1;
                break
            else
                centersS(j,:)=[];
                radiiS(j)=[];
                metricS(j)=[];
                NS=NS-1;
            end
         end
        j=j+1;
    end
    n=n+1;
end



%Remove small particles detections that are inside large particles
NS=length(centersS);
NL=length(centersL);
n=1;
while n<=NS
    j=1;
    while j<=NL
        if sqrt(((centersS(n,1)-centersL(j,1))^2)+((centersS(n,2)-centersL(j,2))^2))<(RS+RL-15)
            if metricS(n)<1.5*metricL(j)
                centersS(n,:)=[];
                radiiS(n)=[];
                metricS(n)=[];
                NS=NS-1;
                n=n-1;
                break
            else
                centersL(j,:)=[];
                radiiL(j)=[];
                metricL(j)=[];
                NL=NL-1;
            end
        end
        j=j+1;
    end
    n=n+1;
end



figure(1)
imshow(img)
viscircles(centersL,radiiL);
viscircles(centersS,radiiS,'Color','b');


