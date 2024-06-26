% Particle Detector, Neighbour finder and Contact validator that creates input files for my Photoelastic Disk Solver
% Particle Detection and Neigbour Finding Adapted from my Earlier Script (joCentersNonMonodisperse.m as of 2016/05/03)
% Photoelastic Disk Solver inspired from peDiskSolve by James Puckett (Phd-Thesis 2012) http://nile.physics.ncsu.edu

% If you use this please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (201X)
% DOI: XXXXXX

% last edit on 2018/08/09 by Joshua Miller (jsmille9@ncsu.edu)

close all % Housekeeping
%clear all % Housekeeping

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           User defined values                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = false; %Generates lots of plots showing results
saveResults = true;

directory = 'C:\Users\Squishfolk\Desktop\Alec\211MSDCF\';
preProDirectory = "C:\Users\Squishfolk\Desktop\Alec\211MSDCF\PreProImages\";
maxFrames = load(fullfile(directory,'maxFrames.mat')); % load maxFrames array from findMaxFrames.m
allFiles = dir(fullfile(directory, 'DSC*.JPG'));
nFrames = length(maxFrames); %how many files are we processing ?

% Hough Transform Values
doParticleDetectionH = true; %Detect particles using Hough Transform?

NsmallH = 256; %Number of small discs
NlargeH = 259; %Number of large discs

% Neighbour Finding Values
findNeighbours = true;

fsigma = 100; %photoelastic stress coefficient
g2cal = 100; %Calibration Value for the g^2 method, can be computed by joG2cal.m

contactG2Threshold = 0.1; %sum of g2 in a contact area larger than this determines a valid contact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 User Input Not Required Below This Line                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use maxFrames to create a list of files to process
% (we only want to process the maxima for now)
files=[];
for n=1:length(allFiles)
    for i=1:length(maxFrames.maxFrames)
        if allFiles(n).name == maxFrames.maxFrames(i)
            filePath = dir(fullfile(directory,allFiles(n).name));
            files=[files; filePath];
        end
    end
end
nFrames=length(files);

for frame = 1:nFrames %Loops for total number of images

    str=sprintf("starting particle detection on frame %d", frame);
    disp(str);


    imageFile = fullfile(directory,files(frame).name); %input filename
    imageFilePrePro = fullfile(preProDirectory,files(frame).name);
    img = imread(imageFile);
    img_prepro = imread(imageFilePrePro);
      
    Rimg_prepro = img_prepro(:,:,:);
    Rimg_prepro = Rimg_prepro*0.4;
    Rimg = img(:,:,1);
    Gimg = img(:,:,2);
    
    %Crop Image to just the portion of the drum that contains particles
    %by drawing circle on first frame in dataset
    if frame==1
    imshow(Rimg_prepro)
    drumCrop = drawcircle; 
    drumMask = createMask(drumCrop);
    
    %Draw a circle around a large particle to determine large particle
    %radius in pixels
    imshow(Rimg_prepro);
    d = drawcircle;
    RL = d.Radius;
    RL = round(RL);
    
    %Draw a circle around a large particle to determine small particle
    %radius in pixels
    d = drawcircle;
    RS = d.Radius;
    RS = round(RS);
    
    pxPerMeter = 0.0047625 / RL;
    end
    

    %Crop images down to particle packing (make the toothed ring and 
    %everything outside of the drum black)
    Rimg_prepro = bsxfun(@times, Rimg_prepro, cast(drumMask, class(Rimg_prepro)));
    wallMask = bsxfun(@times, Gimg, cast(drumMask, class(Gimg)));
    wallMask=imadjust(wallMask,[0 0.05]);
    wallMask=imgaussfilt(wallMask,4);
    wallMask=imbinarize(wallMask);
    Rimg_prepro = bsxfun(@times, Rimg_prepro, cast(wallMask, class(Rimg_prepro)));
    
    if (verbose)
        
        figure(1); %Draw the particle Image
        imshow(img);
        
        figure(2); %Draw the Force Image
        imshow(Gimg);
        
    end
    
    if doParticleDetectionH
        
        particle = AlecDiskFind(Rimg_prepro,Rimg,Gimg,wallMask,pxPerMeter,RS,RL,NsmallH,NlargeH,fsigma);
        
    end
    
    N = length(particle);
    
    if(verbose)
        %add some information about the particles to the plots
        figure(1)
        for n=1:N
            viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
            hold on
            %plot(particle(n).x,particle(n).y,'rx'); %Mark particle centers
            %text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
            hold off
        end
        %{
        figure(2)
        for n=1:N
            viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
            hold on
            plot(particle(n).x,particle(n).y,'rx'); %Mark particle centers
            text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
        end
        %}
        %drawnow;
    end

    % Dim the background of particles in the green channel for proper force
    % fitting
    Gimg = im2double(Gimg);
    Rimg = im2double(Rimg);
    Gimg = Gimg-0.7*Rimg;
    Gimg = Gimg.*(Gimg > 0);
    Gimg = imadjust(Gimg,stretchlim(Gimg));
    Gimg = Gimg - 0.7;
    Gimg = Gimg.*(Gimg > 0);
    
    if findNeighbours

        for n=1:N
        %create a circular mask
        % => Find a better way yo do this masking!
        r = particle(n).r;
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r);
        
        %This crops out a particle
        cropXstart = round(particle(n).x-r);
        cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
        cropYstart = round(particle(n).y-r);
        cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;
        cimg = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImg = cimg.*mask1;     

        particle(n).forceImage=particleImg;
        
        %create a circular mask with a radius that is one pixel smaller
        %for cropping out the relevant gradient

        mask2 = double(sqrt(mask) <= r-1);
        
        %Compute G^2 for each particle
        [gx,gy] = gradient(particleImg);
        g2 = (gx.^2 + gy.^2).*mask2;
        particle(n).g2 = sum(sum(g2));
        particle(n).f = particle(n).g2/g2cal;
        end
        
        dtol = 0.1*RL; % How far away can the outlines of 2 particles be to still be considered Neighbours
        CR = round(0.1*RL); %radius around a contact point that is checked for contact validation
        particle = AlecNeighbourFind(Gimg,wallMask, contactG2Threshold, dtol, CR, verbose, particle);
        
    end
    
    if findNeighbours && verbose
        figure(3);
        imshow(Gimg); 
        hold on
        for n = 1:N
            z = particle(n).z; %get particle coordination number
            if (z>0) %if the particle does have contacts
                for m = 1:z %for each contact
                    %draw contact lines
                    lineX(1)=particle(n).x;
                    lineY(1)=particle(n).y;
                    lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                    lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
                    cX = lineX(1) + (particle(n).r-CR) * cos(particle(n).betas(m));
                    cY = lineY(1) + (particle(n).r-CR) * sin(particle(n).betas(m));
                    hold on; % Don't blow away the image.
                    plot(lineX, lineY,'-y','LineWidth',2);hold on;
                end
            end
        end
    end
    
    if saveResults
    %Save what we got so far
    savePath = fullfile(directory, 'preprocessing files', files(frame).name(1:end-4) + "_preprocessing.mat");
    % Save the file
    save(savePath, 'particle');
    end

   
end