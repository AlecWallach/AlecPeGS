function [particle] = AlecDiskFind3(Rimg_prepro,Gimg,pxPerMeter,RS,RL,NsmallH,NlargeH)


    %Detect large and small particles with generous sensitivity
    [centersL,radiiL,metricL] = imfindcircles(Rimg_prepro,[RL-2 RL+2],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.9775,EdgeThreshold=0.025);
    [centersS,radiiS,metricS] = imfindcircles(Rimg_prepro,[RS-2 RS+2],ObjectPolarity='bright',Method='TwoStage',Sensitivity=0.9775,EdgeThreshold=0.025);
    

    %Create a log of all detections
    Nlarge=length(centersL);
    particleL(1:Nlarge) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[],'meanMin',0,'ringPxVal',0,'meanPxVal',0,'metric',0,'score',0);

    for n=1:Nlarge
        particleL(n).id= n;
        particleL(n).x = centersL(n,1);
        particleL(n).y = centersL(n,2);
        particleL(n).r = radiiL(n);
        particleL(n).metric=metricL(n);
    end
    
    Nsmall=length(centersS);
    particleS(1:Nsmall) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[],'meanMin',0,'ringPxVal',0,'meanPxVal',0,'metric',0,'score',0);

    for n=1:Nsmall
        particleS(n).id= n;
        particleS(n).x = centersS(n,1);
        particleS(n).y = centersS(n,2);
        particleS(n).r = radiiS(n);
        particleS(n).metric=metricS(n);
    end


    %Remove weaker of two overlapping small particles
    NS = length(particleS);
    n=1;
    while n<=NS
        x1=particleS(n).x;
        y1=particleS(n).y;
        m1=particleS(n).metric;
        j=n+1;
        while j<=NS
            x2=particleS(j).x;
            y2=particleS(j).y;
            m2=particleS(j).metric;
            if (sqrt((x1-x2)^2+(y1-y2)^2))<((2*RS)-(.25*RS))
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

    %If two large particlces are overlapping by 15 percent of large radius,
    %remove the one that is less circular (lower metric)
    NL = length(particleL);

    disp(NL)

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
            if sqrt((x1-x2)^2+(y1-y2)^2)<((2*RL)-(.20*RL))
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

    disp(NL)
    
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
            if sqrt((x1-x2)^2+(y1-y2)^2)<(RS+RL-(0.25*RL))
                if m1>3*m2
                    particleL(j)=[];
                    NL=NL-1;
                    j=j-1;
                elseif m1<m2
                    particleS(n)=[];
                    NS=NS-1;
                    n=n-1;
                    break
                end
            end
            j=j+1;
        end
        n=n+1;
    end

    
    disp(NL)

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


    disp(NL)

    %Calculate some values associated with each detection
    %(mean pixel value, minimum pixel value, avg pixel value on outer ring, g2)
    n=1;
    while n<=NL
        r=particleL(n).r;
        x=particleL(n).x;
        y=particleL(n).y;
    
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r);
        mask2 = double(sqrt(mask) <= r-(0.05*r));
        mask3 = double(sqrt(mask) >= r);
        mask4 = mask1.*(double(sqrt(mask) >= r-(0.1*r)));
    
        cropXstart = round(x-r);
        cropXstop = round(x-r)+ size(mask1,1)-1;
        cropYstart = round(y-r);
        cropYstop = round(y-r)+ size(mask1,2)-1;
        cimgG = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
        cimgR = im2double(Rimg_prepro(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImgG = cimgG.*mask1;
        particleImgR = cimgR.*mask1;


        meanPxVal=mean(mean(particleImgG));
        particleL(n).meanPxVal=meanPxVal;

        particleImgMeanMin=particleImgG+mask3;
        meanMin=mean(min(particleImgMeanMin));
        particleL(n).meanMin=meanMin;

        [gx,gy] = gradient(particleImgG);
        g2 = (gx.^2 + gy.^2).*mask2;
        particleL(n).g2 = sum(sum(g2));

        particleImgRing=particleImgR.*mask4;
        particelImgRingMean=mean(mean(particleImgRing));
        particleL(n).ringPxVal=particelImgRingMean;
     
        n=n+1;
    end
    n=1;
    while n<=NS
        r=particleS(n).r;
        x=particleS(n).x;
        y=particleS(n).y;
    
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r);
        mask2 = double(sqrt(mask) <= r-(0.05*r));
        mask3 = double(sqrt(mask) >= r-(0.05*r));
        mask4 = mask1.*(double(sqrt(mask) >= r-(0.1*r)));
    
        cropXstart = round(x-r);
        cropXstop = round(x-r)+ size(mask1,1)-1;
        cropYstart = round(y-r);
        cropYstop = round(y-r)+ size(mask1,2)-1;
        cimgG = im2double(Gimg(cropYstart:cropYstop, cropXstart:cropXstop));
        cimgR = im2double(Rimg_prepro(cropYstart:cropYstop, cropXstart:cropXstop));
        particleImgG = cimgG.*mask1;
        particleImgR = cimgR.*mask1;

        particleImgMeanMin=particleImgG+mask3;
        meanMin=mean(min(particleImgMeanMin));
        particleS(n).meanMin=meanMin;

        [gx,gy] = gradient(particleImgG);
        g2 = (gx.^2 + gy.^2).*mask2;
        particleS(n).g2 = sum(sum(g2));

        particleImgRing=particleImgR.*mask4;
        particelImgRingMean=mean(mean(particleImgRing));
        particleS(n).ringPxVal=particelImgRingMean;

        n=n+1;
    end
    

    %If average of minimum pixel value in each col of particle is very low,
    %it is on the toothed ring so we delete it
    n=1;
    while n<=NS
        if particleS(n).meanMin < 0.075
            particleS(n) = [];
            NS=NS-1;
            n=n-1;
        end
        n=n+1;
    end
    n=1;
    while n<=NL
        if particleL(n).meanMin < 0.05
            particleS(n) = [];
            NL=NL-1;
            n=n-1;
        end
        n=n+1;
    end

    ringPxValsL = zeros(NL,1);
    for n = 1:NL
        ringPxValsL = (particleL(n).ringPxVal);
    end
    while NL > NlargeH
        [~,nMin] = min(ringPxValsL);
        particleL(nMin)=[];
        ringPxValsL(nMin) = [];
        NL = NL - 1;
    end

    ringPxValsS = zeros(NS,1);
    for n = 1:NS
        ringPxValsS = (particleS(n).ringPxVal);
    end
    while NS > NsmallH
        [~,nMin] = min(ringPxValsS);
        particleS(nMin)=[];
        ringPxValsS(nMin) = [];
        NS = NS - 1;
    end

    

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


    particle=[particleL, particleS];
    N=length(particle);
    for n=1:N
        particle(n).id = n;
        particle(n).rm = particle(n).r * pxPerMeter;
        if particle(n).r < ( mean(RL)+ mean(RS) ) /2
            particle(n).color = 'b';
        else
            particle(n).color = 'r';
        end
    end

end


