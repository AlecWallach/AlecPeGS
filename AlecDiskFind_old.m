function [particle,meanListL] = AlecDiskFind_old(Rimg,Rimg_prepro,Gimg,pxPerMeter,fsigma,RS,RL,NsmallH,NlargeH)


    %Detect large and small particles with generous sensitivity
    [centersL,radiiL,metricL] = imfindcircles(Rimg_prepro,[RL-2 RL],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.9875,EdgeThreshold=0.0225);
    [centersS,radiiS,metricS] = imfindcircles(Rimg_prepro,[RS-2 RS],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.9775,EdgeThreshold=0.025);
    
    
    %Remove weaker of two overlapping small particles
    NS = length(centersS);
    n=1;
    while n<=NS
        j=n+1;
        while j<=NS
            if (sqrt(((centersS(n,1)-centersS(j,1))^2)+((centersS(n,2)-centersS(j,2))^2))<((2*RS)-(.25*RS)))
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
                    j=j-1;
                end
             end
            j=j+1;
        end
        n=n+1;
    end
    
    
    %If two large particlces are overlapping by 15 percent of large radius,
    %remove the one that is less circular (lower metric)
    NL = length(centersL);
    n=1;
    while n<=NL
        j=n+1;
        while j<=NL
            if (sqrt(((centersL(n,1)-centersL(j,1))^2)+((centersL(n,2)-centersL(j,2))^2))<((2*RL)-(.15*RL)))
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
                    j=j-1;
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
            if sqrt(((centersS(n,1)-centersL(j,1))^2)+((centersS(n,2)-centersL(j,2))^2))<(RS+RL-(0.25*RS))
                if metricS(n)<2.5*metricL(j)
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
                    j=j-1;
                end
            end
            j=j+1;
        end
        n=n+1;
    end

 
    %If the detected particle doesn't have many low value green pixels, its
    %probably on the teeth, so remove it
    %{
    NS=length(centersS);
    meanListS=[];
    n=1;
    while n<=NS
        r=radiiS(n);
        x=centersS(n,1);
        y=centersS(n,2);
    
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r-0.05*r);
        mask2 = double(sqrt(mask) >= r-0.05*r);
    
        cropXstart = round(x-r);
        cropXstop = round(x-r)+ size(mask1,1)-1;
        cropYstart = round(y-r);
        cropYstop = round(y-r)+ size(mask1,2)-1;
        cimg = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImg = cimg.*mask1;
        particleImg=particleImg+mask2;

        meanMin=mean(min(particleImg));
        meanListS=[meanListS meanMin];

        if meanMin>0.6
            centersS(n,:)=[];
            radiiS(n)=[];
            metricS(n)=[];
            NS=NS-1;
            n=n-1;
        end
        n=n+1;
    end
    %}
    NL=length(centersL);
    n=1;
    meanListL=[];
    while n<=NL
        r=radiiL(n);
        x=centersL(n,1);
        y=centersL(n,2);
    
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r-0.05*r);
        mask2 = double(sqrt(mask) >= r-0.05*r);
    
        cropXstart = round(x-r);
        cropXstop = round(x-r)+ size(mask1,1)-1;
        cropYstart = round(y-r);
        cropYstop = round(y-r)+ size(mask1,2)-1;
        cimg = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImg = cimg.*mask1;
        particleImg=particleImg+mask2;
    
        meanMin=mean(min(particleImg));
        meanListL=[meanListL meanMin];
        
        if meanMin>0.47
            centersL(n,:)=[];
            radiiL(n)=[];
            metricL(n)=[];
            NL=NL-1;
            n=n-1;
        end
        
        n=n+1;
    end
    
    
    %If a particle has less than 1 neighbor, remove it
    NS=length(centersS);
    NL=length(centersL);
    n=1;
    while n<=NS
        j=1;
        z=0;
        while j<=NS
            if sqrt(((centersS(n,1)-centersS(j,1))^2)+((centersS(n,2)-centersS(j,2))^2))<((2*RS)+(0.20*RS)) && n~=j
                z=z+1;
            end
            j=j+1;
        end
        i=1;
        while i<=NL
            if sqrt(((centersS(n,1)-centersL(i,1))^2)+((centersS(n,2)-centersL(i,2))^2))<((RS+RL)+(0.20*RL))
                z=z+1;
            end
            i=i+1;
        end
        if z<1
            centersS(n,:)=[];
            radiiS(n)=[];
            metricS(n)=[];
            NS=NS-1;
            n=n-1;
        end
        n=n+1;
    end
    NS=length(centersS);
    NL=length(centersL);
    n=1;
    while n<=NL
        j=1;
        z=0;
        while j<=NL
            if sqrt(((centersL(n,1)-centersL(j,1))^2)+((centersL(n,2)-centersL(j,2))^2))<((2*RL)+(0.20*RL)) && n~=j
                z=z+1;
            end
            j=j+1;
        end
        i=1;
        while i<=NS
            if sqrt(((centersL(n,1)-centersS(i,1))^2)+((centersL(n,2)-centersS(i,2))^2))<((RS+RL)+(0.20*RL))
                z=z+1;
            end
            i=i+1;
        end
        if z<1
            centersL(n,:)=[];
            radiiL(n)=[];
            metricL(n)=[];
            NL=NL-1;
            n=n-1;
        end
        n=n+1;
    end
    
    
    %If the mean pixel value inside a particle is very low, remove it 
    %(it's likely a detection outside of the drum)
    NS=length(centersS);
    n=1;
    particleImgMeanS=[];
    while n<=NS
        r=radiiS(n);
        x=centersS(n,1);
        y=centersS(n,2);
    
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r);
        mask2 = double(sqrt(mask) >= r-(0.1*RS));
        mask = mask1.*mask2;
    
        cropXstart = round(x-r);
        cropXstop = round(x-r)+ size(mask,1)-1;
        cropYstart = round(y-r);
        cropYstop = round(y-r)+ size(mask,2)-1;
        cimg = im2double(Rimg(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImg = cimg.*mask;
    
        meanPxVal=mean(mean(particleImg));
        particleImgMeanS = [particleImgMeanS; meanPxVal];
    
        n=n+1;
    end
    while NS>NsmallH
        [~,nMin]=min(particleImgMeanS);
        centersS(nMin,:)=[];
        radiiS(nMin)=[];
        metricS(nMin)=[];
        particleImgMeanS(nMin)=[];
        NS=NS-1;
    end
    
    NL=length(centersL);
    n=1;
    particleImgMeanL=[];
    while n<=NL
        r=radiiL(n);
        x=centersL(n,1);
        y=centersL(n,2);
    
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r);
        mask2 = double(sqrt(mask) >= r-(0.1*RL));
        mask = mask1.*mask2;
    
        cropXstart = round(x-r);
        cropXstop = round(x-r)+ size(mask,1)-1;
        cropYstart = round(y-r);
        cropYstop = round(y-r)+ size(mask,2)-1;
        cimg = im2double(Rimg(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImg = cimg.*mask;
    
        meanPxVal=mean(mean(particleImg));
        particleImgMeanL = [particleImgMeanL; meanPxVal];
     
        n=n+1;
    end
    while NL>NlargeH
        [~,nMin]=min(particleImgMeanL);
        centersL(nMin,:)=[];
        radiiL(nMin)=[];
        metricL(nMin)=[];
        particleImgMeanL(nMin)=[];
        NL=NL-1;
    end
    
    
    %Bookkeeping
    Nlarge=length(centersL);
    particle(1:Nlarge) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
    
    for n=1:Nlarge
        particle(n).id= n;
        particle(n).x = centersL(n,1);
        particle(n).y = centersL(n,2);
        particle(n).r = radiiL(n);
    end
    
    Nsmall=length(centersS);
    N = Nsmall + Nlarge;
    particle(Nlarge+1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
    
    %Bookkeeping
    for n=1:Nsmall
        particle(Nlarge+n).id= Nlarge+n;
        particle(Nlarge+n).x = centersS(n,1);
        particle(Nlarge+n).y = centersS(n,2);
        particle(Nlarge+n).r = radiiS(n);
    end
    
    for n=1:N
        particle(n).rm = particle(n).r * pxPerMeter;
        particle(n).fsigma = fsigma;
        if particle(n).r < ( mean(RL)+ mean(RS) ) /2
            particle(n).color = 'b';
        else
            particle(n).color = 'r';
        end
    end

end
