function [xyz, front_area, weight] = gen_payload()
% 100 g bags only
vol = 0;
while vol < 3 
    xlim = 3;
    ylim = 1;
    zlim = 9;
    x = randi([1 xlim],1);
    y = randi([1,ylim],1);
    z = randi([1,zlim],1);
    vol = x*y*z;
end
xyz = [x y z];
front_area = y*z*0.00165; %m^2
if front_area < 0.00894
    front_area = 0.00894;
end
if front_area > 0.019081
    front_area = 0.019081;
end

weight = 9.81*0.1*vol;

end