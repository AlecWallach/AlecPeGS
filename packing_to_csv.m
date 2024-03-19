input_file = load('/Volumes/SD/DCIM/212MSDCF/packingStruct.mat');

packing_struct = input_file.packing;

angle = zeros(60,1);
for n=1:60
    angle(n) = packing_struct(n).packingAngleDegrees;
end
time = 1:60;

plot(time, angle, LineWidth=2)
xlabel("Time (s)")
ylabel("Angle (degrees)")
title("Packing Angle vs. Time")
