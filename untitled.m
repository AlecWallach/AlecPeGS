z=2; %Number of contacts this particle has
%f = [0.6 0.6]; %Absolute forces on this particle
alpha = [0 0];  %Alpha contact angles on this particle
beta = [1.7329 0.8144];  %Beta contact angles on this particle
fsigma = 65;  %Photoelastic stress coefficient of this particle
rm = 0.0047625; %Particle radius in meters
px = 150; %Return image size in pixels
f = [0.0192 0.0164];

img = joForceImg (z, f, alpha, beta, fsigma, rm, px, verbose);

imshow(img)


