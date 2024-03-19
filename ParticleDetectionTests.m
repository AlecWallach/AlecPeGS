img = imread('/Users/alecwallach/Desktop/Squishlab/Imaging/DSC_3676.JPG');
Rimg = img(:,:,1);

figure(1)
imshow(Rimg);
particle = AlecDiskFind(Rimg,0.98,0.1);

N = length(particle);

figure(1)
for n=1:N
    viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color); %draw particle outline
    hold on
    plot(particle(n).x,particle(n).y,'rx'); %Mark particle centers
    text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','w');
end