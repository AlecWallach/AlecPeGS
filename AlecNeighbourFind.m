function [particle] = AlecNeighbourFind(Gimg, wallMask, contactG2Threshold, dtol, CR, verbose, particle)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

N = length(particle);

xmat = zeros([N,1]);
ymat = zeros([N,1]); % Preallocation
rmat = zeros([N,1]);

for l = 1:N
    xmat(l) = particle(l).x;
    ymat(l) = particle(l).y; %Pulls data from particle structure
    rmat(l) = particle(l).r;
end

rmats = rmat; %Saves our radius matrix for later

dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
rmat = rmat + rmat'; %Makes a combination of radii for each particle

friendmat = dmat < (rmat + dtol) & dmat~=0; %Logical "friend" matrix

friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
[f1, f2] = find(friendmat == 1); %Creates an index of particles that are considered touching

%Create circular mask with radius equal to CR
maskCR = abs(-CR:CR);
maskCR = maskCR.^2 + maskCR.^2';
maskCR = double(sqrt(maskCR) <= CR-1);

for l = 1:length(f1)
    x1 = particle(f1(l)).x;
    y1 = particle(f1(l)).y;
    r1 = particle(f1(l)).r;
    x2 = particle(f2(l)).x;
    y2 = particle(f2(l)).y;
    r2 = particle(f2(l)).r;
    
    contactXp1 = x1 + (r1 - CR) * cos(atan2(y2-y1,x2-x1));
    contactYp1 = y1 + (r1 - CR) * sin(atan2(y2-y1,x2-x1));
    
    contactXp2 = x1 + (r1 + CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * cos(atan2(y2-y1,x2-x1));
    contactYp2 = y1 + (r1 + CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * sin(atan2(y2-y1,x2-x1));
    
    contactImg = im2double(imcrop(Gimg,[contactXp1-CR contactYp1-CR CR*2 CR*2]));
    contactImg = contactImg.*maskCR;
    
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2);
    contactG2p1 = sum(sum(g2));
    contactIp1 = sum(sum(contactImg));
    
    contactImg = im2double(imcrop(Gimg,[contactXp2-CR contactYp2-CR CR*2 CR*2]));
    contactImg = contactImg.*maskCR;
    
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2);
    contactG2p2 = sum(sum(g2));
    contactIp2 = sum(sum(contactImg));
    
    %if we declare our contact valid
    %if(contactG2p1 > contactG2Threshold && contactG2p2 > contactG2Threshold)
    %cI = sum(sum(contactImg)); %Use integrated intensity instead of g2
    %if(cI > contactIThreshold) %Use integrated intensity instead of g2
    %Plot contact area
    if (verbose)
        display(['contact found between particles ',num2str(f1(l)),' and ',num2str(f2(l))]);
        %viscircles([contactXp2; contactYp2]', CR,'EdgeColor','w');
        %text(contactXp1, contactYp1,num2str(contactG2p1),'Color','w');
        %drawnow;
    end
    %this is a valid contact, remember it
    particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
    particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p1; %remember the g2 value of the current contact area
    particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp1;
    particle(f1(l)).neighbours(particle(f1(l)).z) = f2(l); %particle m is now noted as a neigbour in the particle l datastructure
    particle(f1(l)).betas(particle(f1(l)).z) = atan2(y2-y1,x2-x1); %the contact angle to particle m is now noted in the particle l datastructure
    particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
    particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p2; %remember the g2 value of the current contact area
    particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp2;
    particle(f2(l)).neighbours(particle(f2(l)).z) = f1(l); %particle m is now noted as a neigbour in the particle l datastructure
    particle(f2(l)).betas(particle(f2(l)).z) = atan2(y1-y2,x1-x2);
    %end
end


%Check for contacts with the wall
circs = [xmat,ymat,rmats];
for n=1:length(circs)
    x=circs(n,1);
    y=circs(n,2);
    r=circs(n,3);
    r=r+0.1*r;

    %Find all of the points on the perimeter of the particle near a
    %wall
    possibleContactAngles = [];
    for theta = linspace(3*pi/2,(7*pi/2)-(2*pi/100),100)

        x_test = round(x + r * cos(theta));
        y_test = round(y + r * sin(theta)); 

        if wallMask(y_test, x_test) == 0
            possibleContactAngles = [possibleContactAngles; theta];
        end
    end

    if ~isempty(possibleContactAngles)

        %Determine if the particle has 1 or 2 contacts with the wall
        %(in a corner or not in a corner)

        possibleContactAngles = sort(possibleContactAngles);

        multipleContacts=false;
        for j=1:length(possibleContactAngles)-1
            if possibleContactAngles(j+1)-possibleContactAngles(j) > pi/5
                multipleContacts = true;
                break
            end
        end
    
        %The recorded contact angle is the average angle of the possible
        %contacts for a continuous section
        if multipleContacts
            contactAngles = [mean(possibleContactAngles(1:j)); mean(possibleContactAngles(j+1:end))];
        else
            contactAngles = (mean(possibleContactAngles));
        end
    
        %Record contact, contact angle, g2
        for k=1:length(contactAngles)
            contactX = round(x + r * cos(contactAngles(k)));
            contactY = round(y + r * cos(contactAngles(k)));
    
            contactImg = im2double(imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
            contactImg = contactImg.*maskCR;
    
            [gx,gy] = gradient(contactImg);
            g2 = (gx.^2 + gy.^2);
            contactG2 = sum(sum(g2));
            
            
            cI = sum(sum(contactImg));
            %if(cI > contactIThreshold)
            %this is a valid contact, remember it
            if(verbose)
                %text(contactX,contactY,num2str(contactG2),'Color','w');
                %viscircles([contactX; contactY]', CR,'EdgeColor','w');
            end
            particle(n).z= particle(n).z +1; %increase coordination number
            particle(n).contactG2s(particle(n).z)=contactG2;
            particle(n).contactIs(particle(n).z)=cI;
            particle(n).neighbours(particle(n).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
            particle(n).betas(particle(n).z) = contactAngles(k); %the contact angle to the wall is now noted in the particle l datastructure
        end
    end
        
end