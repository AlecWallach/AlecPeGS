%test
function [particle] = AlecDiskFind(Rimg_prepro,Rimg,Gimg,wallMask,pxPerMeter,RS,RL,NsmallH,NlargeH,fsigma)


    %Detect large and small particles with generous sensitivity
    [centersL,radiiL,metricL] = imfindcircles(Rimg_prepro,[RL-3 RL+3],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.97,EdgeThreshold=0.025);
    [centersS,radiiS,metricS] = imfindcircles(Rimg_prepro,[RS-2 RS+2],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.97,EdgeThreshold=0.025);
    
    fprintf('imfindcircles initially detected %d large particles\n',length(centersL))
    fprintf('imfindcircles initially detected %d small particles\n',length(centersS))

    %Create a log of all detections
    Nlarge=length(centersL);
    particleL(1:Nlarge) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[],'meanMin',0,'ringPxVal',0,'meanPxVal',0,'metric',0,'ptclIntensity',0,'score',0);

    for n=1:Nlarge
        particleL(n).id= n;
        particleL(n).x = centersL(n,1);
        particleL(n).y = centersL(n,2);
        particleL(n).r = radiiL(n);
        particleL(n).metric=metricL(n);
    end
    
    Nsmall=length(centersS);
    particleS(1:Nsmall) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[],'meanMin',0,'ringPxVal',0,'meanPxVal',0,'metric',0,'ptclIntensity',0,'score',0);

    for n=1:Nsmall
        particleS(n).id= n;
        particleS(n).x = centersS(n,1);
        particleS(n).y = centersS(n,2);
        particleS(n).r = radiiS(n);
        particleS(n).metric=metricS(n);
    end

    %Remove any small or large particle whose perimeter lies on the
    %toothed ring
    NS=length(particleS);
    n = 1;
    while n <= NS
        x = particleS(n).x;
        y = particleS(n).y;
        r = 0.67 * particleS(n).r;
 
    
        for theta = linspace(0,2*pi)
            x_test = round(x + r * cos(theta));
            y_test = round(y + r * sin(theta));
    
            if wallMask(y_test, x_test) == 0
                particleS(n)=[];
                NS=NS-1;
                n=n-1;
                break
            end
        end
        n=n+1;
    end
    NL=length(particleL);
    n=1;
    while n<=NL
        x=particleL(n).x;
        y=particleL(n).y;
        r=0.67*(particleL(n).r);

        for theta = linspace(0,2*pi)
            x_test=round(x+r*cos(theta));
            y_test=round(y+r*sin(theta));

            if wallMask(y_test,x_test)==0
                particleL(n)=[];
                NL=NL-1;
                n=n-1;
                break;
            end
        end
        n=n+1;
    end

    fprintf('%d large particles remaining after removing particles on boundary\n',NL)
    fprintf('%d small particles remaining after removing particles on boundary\n',NS)

    
    %If two small particles are overlapping by more than 20% of a small
    %particle radius, remove the particle with lower metric
    n=1;
    while n<=NS
        x1=particleS(n).x;
        y1=particleS(n).y;
        r1=particleS(n).r;
        m1=particleS(n).metric;
        j=n+1;
        while j<=NS
            x2=particleS(j).x;
            y2=particleS(j).y;
            r2=particleS(j).r;
            m2=particleS(j).metric;
            if (sqrt((x1-x2)^2+(y1-y2)^2))<(0.5*r1+r2)
                if(m1<m2)
                    particleS(n)=[];
                    NS=NS-1;
                    n=n-1;
                    break
                else
                    particleS(j)=[];
                    NS=NS-1;
                    j=j-1;
                end
             end
            j=j+1;
        end
        n=n+1;
    end

    fprintf('%d small particles remaining after removing overlapping small particles\n',NS)

    %If two large particlces are overlapping by 20% of large radius,
    %remove the one that is less circular (lower metric)
    NL = length(particleL);
    n=1;
    while n<=NL
        x1=particleL(n).x;
        y1=particleL(n).y;
        m1=particleL(n).metric;
        j=n+1;
        while j<=NL
            x2=particleL(j).x;
            y2=particleL(j).y;
            m2=particleL(j).metric;
            if sqrt((x1-x2)^2+(y1-y2)^2)<((2*RL)-(.3*RL))
                if(m1<m2)
                    particleL(n)=[];
                    NL=NL-1;
                    n=n-1;
                    break
                else
                    particleL(j)=[];
                    NL=NL-1;
                    j=j-1;
                end
             end
            j=j+1;
        end
        n=n+1;
    end

    fprintf('%d large particles remaining after removing overlapping large particles\n',NL)

    %If a large particle is enclosing multiple small particles at this
    %point it can't be a real particle
    n=1;
    while n<=NL
        x1=particleL(n).x;
        y1=particleL(n).y;
        r1=particleL(n).r;
        m1=particleL(n).metric;
        j=1;
        numParticlesPartialEnclosed=0;
        numParticlesFullyEnclosed=0;
        while j<=NS
            x2=particleS(j).x;
            y2=particleS(j).y;
            r2=particleS(j).r;
            m2=particleS(j).metric;
            if sqrt((x1-x2)^2+(y1-y2)^2) < r1+0.5*r2 && m2>1.25*m1
                numParticlesPartialEnclosed=numParticlesPartialEnclosed+1;
            end
            if sqrt((x1-x2)^2+(y1-y2)^2) < r1-0.8*r2 && m2>1.25*m1
                numParticlesFullyEnclosed=numParticlesFullyEnclosed+1;
            end
            j=j+1;
        end
        if numParticlesPartialEnclosed>1
            particleL(n)=[];
            NL=NL-1;
            n=n-1;
        elseif numParticlesFullyEnclosed>0
            particleL(n)=[];
            NL=NL-1;
            n=n-1;
        end
        n=n+1;
    end

    
    %Remove small particles detections that are inside large particles
    n=1;
    while n<=NS
        x1=particleS(n).x;
        y1=particleS(n).y;
        m1=particleS(n).metric;
        j=1;
        while j<=NL
            x2=particleL(j).x;
            y2=particleL(j).y;
            m2=particleL(j).metric;
            if sqrt((x1-x2)^2+(y1-y2)^2)<RL-0.75*RS && m1<2*m2
                    particleS(n)=[];
                    NS=NS-1;
                    n=n-1;
                    break
            end
            j=j+1;
        end
        n=n+1;
    end
    

    fprintf('%d large particles remaining after removing large particles overlapping small\n',NL)
    fprintf('%d small particles remaining after removing small particles inside large\n',NS);


    %Remove small particle false positive detections that are sitting in
    %empty voids
    n=1;
    while n<=NS
        x1=particleS(n).x;
        y1=particleS(n).y;
        m1=particleS(n).metric;
        j=1;
        numLargeOverlap=0;
        numSmallOverlap=0;
        while j<=NL
            x2=particleL(j).x;
            y2=particleL(j).y;
            m2=particleL(j).metric;
            if sqrt((x1-x2)^2+(y1-y2)^2)<(RL+0.5*RS) && m1<m2
                numLargeOverlap=numLargeOverlap+1;
            end
            if sqrt((x1-x2)^2+(y1-y2)^2)<(RL+0.8*RS) && m1<m2
                numSmallOverlap=numSmallOverlap+1;
            end
            if numLargeOverlap>0 || numLargeOverlap+numSmallOverlap>1
                particleS(n)=[];
                    NS=NS-1;
                    n=n-1;
                    break
            end
            j=j+1;
        end
        %{
        k=1;
        while k<=NS
            x3=particleS(k).x;
            y3=particleS(k).y;
            m3=particleS(k).metric;
            if sqrt((x1-x3)^2+(y1-y3)^2)<(RS+0.6*RS) && m1<m3
                numSmallOverlap=numSmallOverlap+1;
            end
            if numLargeOverlap>0 || numLargeOverlap+numSmallOverlap>1
                particleS(n)=[];
                    NS=NS-1;
                    n=n-1;
                    break
            end
            k=k+1;
        end
        %}
        n=n+1;
    end

    fprintf('%d small particles remaining after removing small particles overlapping large\n',NS)

    %If a particle has less than 1 neighbor, remove it
    n=1;
    while n<=NS
        x1=particleS(n).x;
        y1=particleS(n).y;
        id1=particleS(n).id;
        j=1;
        z=0;
        while j<=NS
            x2=particleS(j).x;
            y2=particleS(j).y;
            id2=particleS(j).id;
            if sqrt((x1-x2)^2+(y1-y2)^2)<((2*RS)+(0.20*RS)) && id1~=id2
                z=z+1;
            end
            j=j+1;
        end
        i=1;
        while i<=NL
            x2=particleL(i).x;
            y2=particleL(i).y;
            if sqrt((x1-x2)^2+(y1-y2)^2)<((RS+RL)+(0.20*RL))
                z=z+1;
            end
            i=i+1;
        end
        if z<1
            particleS(n)=[];
            NS=NS-1;
            n=n-1;
        end
        n=n+1;
    end
    n=1;
    while n<=NL
        x1=particleL(n).x;
        y1=particleL(n).y;
        id1=particleL(n).id;
        j=1;
        z=0;
        while j<=NL
            x2=particleL(j).x;
            y2=particleL(j).y;
            id2=particleL(j).id;
            if sqrt((x1-x2)^2+(y1-y2)^2)<((2*RL)+(0.20*RL)) && id1~=id2
                z=z+1;
            end
            j=j+1;
        end
        i=1;
        while i<=NS
            x2=particleS(i).x;
            y2=particleS(i).y;
            if sqrt((x1-x2)^2+(y1-y2)^2)<((RS+RL)+(0.20*RL))
                z=z+1;
            end
            i=i+1;
        end
        if z<1
            particleL(n)=[];
            NL=NL-1;
            n=n-1;
        end
        n=n+1;
    end

    fprintf('%d large particles remaining after trimming isolated particles\n',NL)
    fprintf('%d small particles remaining after trimming isolated particles\n',NS)



    %The following section scores the remaining particles based on several
    %criteria and wittles down the detections to the correct number of
    %particles

    if NL>NlargeH
        n=1;
        while n<=NL
            r=particleL(n).r;
            x=particleL(n).x;
            y=particleL(n).y;
        
            mask = abs(-r:r);
            mask = mask.^2 + mask.^2';
            mask1 = double(sqrt(mask) <= r);
            mask2 = double(sqrt(mask) <= r);
            mask3 = double(sqrt(mask) >= r);
            mask4 = mask1.*(double(sqrt(mask) >= r-(0.1*r)));
        
            cropXstart = round(x-r);
            cropXstop = round(x-r)+ size(mask1,1)-1;
            cropYstart = round(y-r);
            cropYstop = round(y-r)+ size(mask1,2)-1;
            cimgG = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
            cimgR = im2double(Rimg(cropYstart:cropYstop, cropXstart:cropXstop));
            particleImgG = cimgG.*mask1;
            particleImgR = cimgR.*mask1;
    
            %mean pixel value for the particle in the red channel
            meanPxVal=mean(mean(particleImgR));
            particleL(n).meanPxVal=meanPxVal;

            %particle intensity in the red channel
            ptclIntensity=sum(sum(particleImgR));
            particleL(n).ptclIntensity=ptclIntensity;
    
            %average value of the minimum pixel in each column of the particle
            %in the green channel
            particleImgMeanMin=particleImgG+mask3;
            meanMin=mean(min(particleImgMeanMin));
            particleL(n).meanMin=meanMin;
    
            %the g2 value of the particle
            [gx,gy] = gradient(particleImgG);
            g2 = (gx.^2 + gy.^2).*mask2;
            particleL(n).g2 = sum(sum(g2));
    
            %the average value of the pixels in the outermost 10% of the
            %particle
            particleImgRing=particleImgR.*mask4;
            particelImgRingMean=mean(mean(particleImgRing));
            particleL(n).ringPxVal=particelImgRingMean;
         
            n=n+1;
        end

        %Find medians of the above quantites
        scoresL=zeros(NL,1);
        ringPxValL=zeros(NL,1);
        G2L=zeros(NL,1);
        ptclIntensityL=zeros(NS,1);
        meanPxValL=zeros(NL,1);
        meanMinL=zeros(NL,1);
        for n=1:NL
            ringPxValL(n)=particleL(n).ringPxVal;
            meanPxValL(n)=particleL(n).meanPxVal;
            ptclIntensityL(n)=particleL(n).ptclIntensity;
            meanMinL(n)=particleL(n).meanMin;
            G2L(n)=particleL(n).g2;
        end
        avgRingPxValL=median(ringPxValL);
        avgMeanPxValL=median(meanPxValL);
        meanPtclIntensityL=mean(ptclIntensityL);
        avgMeanMinL=median(meanMinL);
        avgG2L=median(G2L);

        %Give each particle a score based on deviation from mean values
        for n=1:NL
            meanMin=particleL(n).meanMin;
            g2=particleL(n).g2;
            ringPxVal=particleL(n).ringPxVal;
            ptclIntensity=particleL(n).ptclIntensity;
            metric=particleL(n).metric;
            meanPxVal=particleL(n).meanPxVal;
            
            %Score
            particleL(n).score = abs(ringPxVal-avgRingPxValL)+abs(meanPxVal-avgMeanPxValL);
            scoresL(n)=particleL(n).score;
        end

        %Make all scores positive and sort
        minScoreL=abs(min(scoresL));
        for n=1:NL
            particleL(n).score = particleL(n).score + minScoreL;
            scoresL(n)=particleL(n).score;
        end
        scoresL = sort(scoresL);

        %Remove particles with highest scores until we reach known number of
        %particles
        n=1;
        while n<=NL
            if particleL(n).score > scoresL(NlargeH)            
                particleL(n) = [];
                n=n-1;
                NL=NL-1;
                if length(particleL)==NlargeH
                    break
                end
                
            end
            n=n+1;
        end

    end


    if NS>NsmallH
        n=1;
        while n<=NS
            r=particleS(n).r;
            x=particleS(n).x;
            y=particleS(n).y;
        
            mask = abs(-r:r);
            mask = mask.^2 + mask.^2';
            mask1 = double(sqrt(mask) <= 0.4*r);
            mask2 = double(sqrt(mask) <= r);
            mask3 = double(sqrt(mask) >= r-(0.05*r));
            mask4 = mask1.*(double(sqrt(mask) >= r-(0.1*r)));
        
            cropXstart = round(x-r);
            cropXstop = round(x-r)+ size(mask1,1)-1;
            cropYstart = round(y-r);
            cropYstop = round(y-r)+ size(mask1,2)-1;
            cimgG = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
            cimgR = im2double(Rimg(cropYstart:cropYstop, cropXstart:cropXstop));
            particleImgG = cimgG.*mask1;
            particleImgR = cimgR.*mask1;
    
            %mean pixel value for the particle in the red channel
            meanPxVal=mean(mean(particleImgR));
            particleS(n).meanPxVal=meanPxVal;

            %particle intensity in the red channel
            ptclIntensity=sum(sum(particleImgR));
            particleS(n).ptclIntensity=ptclIntensity;
            
            %average value of the minimum pixel in each column of the particle
            %in the green channel
            particleImgMeanMin=particleImgG+mask3;
            meanMin=mean(min(particleImgMeanMin));
            particleS(n).meanMin=meanMin;
    
            %the g2 value of the particle
            [gx,gy] = gradient(particleImgG);
            g2 = (gx.^2 + gy.^2).*mask2;
            particleS(n).g2 = sum(sum(g2));
    
            %the average value of the pixels in the outermost 10% of the
            %particle
            particleImgRing=particleImgR.*mask4;
            particelImgRingMean=mean(mean(particleImgRing));
            particleS(n).ringPxVal=particelImgRingMean;
    
            n=n+1;
        end

        scoresS=zeros(NS,1);
        ringPxValS=zeros(NS,1);
        meanPxValS=zeros(NS,1);
        ptclIntensityS=zeros(NS,1);
        meanMinS=zeros(NS,1);
        G2S=zeros(NS,1);
        for n=1:NS
            ringPxValS(n)=particleS(n).ringPxVal;
            meanPxValS(n)=particleS(n).meanPxVal;
            ptclIntensityS(n)=particleS(n).ptclIntensity;
            meanMinS(n)=particleS(n).meanMin;
            G2S(n)=particleS(n).g2;
        end
        avgRingPxValS=mean(ringPxValS);
        avgMeanPxValS=mean(meanPxValS);
        meanPtclIntensityS=median(ptclIntensityS);
        avgMeanMinS=mean(meanMinS);
        avgG2S=mean(G2S);

        %Give each particle a score based on deviation from mean values
        for n=1:NS
            meanMin=particleS(n).meanMin;
            g2=particleS(n).g2;
            ringPxVal=particleS(n).ringPxVal;
            ptclIntensity=particleS(n).ptclIntensity;
            metric=particleS(n).metric;
            meanPxVal=particleS(n).meanPxVal;
            
            %Score
            particleS(n).score = -(ptclIntensity-meanPtclIntensityS);
            scoresS(n)=particleS(n).score;
        end
        
        %Make all scores positive and sort
        minScoreS=abs(min(scoresS));
        for n=1:NS
        particleS(n).score = particleS(n).score + minScoreS;
        scoresS(n)=particleS(n).score;
        end
        scoresS = sort(scoresS);


        %Remove particles with highest scores until we reach known number
        %of particles
        n=1;
        while n<=NS
            if particleS(n).score > scoresS(NsmallH)
                particleS(n) = [];
                n=n-1;
                NS=NS-1;
                if length(particleS)==NsmallH
                    break
                end
            end
            n=n+1;
        end
    end
    

    %Combine large and small particles into a single struct
    particle=[particleL, particleS];
    N=length(particle);
    for n=1:N
        particle(n).id = n;
        particle(n).rm = particle(n).r * pxPerMeter;
        particle(n).fsigma = fsigma;
        if particle(n).r < ( mean(RL)+ mean(RS) ) /2
            particle(n).color = 'b';
        else
            particle(n).color = 'r';
        end
    end
end


