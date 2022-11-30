
close all;
clear all;


[X,Y] = meshgrid(4:0.01:5.5,5:0.01:8.5);
Z = sin(X) + cos(Y)+2;
surf(X,Y,Z, "FaceColor", "flat", "FaceAlpha", 0.7, "EdgeColor", "none")
hold on;

[data_X, data_Y] = meshgrid(4.5:0.5:5,6.25:1:7.25);
data_Z = sin(data_X) + cos(data_Y)+2;
surf(data_X,data_Y,data_Z, "FaceColor", "none" );
hold on;
% project corners down to xy plane
for i = 1:length(data_X(1,:))
    for j = 1:length(data_X(:,1))
        plot3([data_X(i,j), data_X(i,j)], ...
              [data_Y(i,j), data_Y(i,j)], ...
              [0,data_Z(i,j)], ...
              '--', "Color", [0.6350 0.0780 0.1840]);
        hold on;
        plot3([data_X(i,j)], ...
              [data_Y(i,j)], ...
              [0], ...
              'o', "Color",  [0.6350 0.0780 0.1840], ...
              "MarkerFaceColor", [0.6350 0.0780 0.1840])
        hold on;
    end
end


% plot sampled point
x = 4.8;
y = 6.7;
z = bilinearly_interpolate(data_X(1,:), data_Y(:,1), data_Z.', x, y);
plot3([x],[y],[z], ...
      'o', "Color",  "#FF8000", ...
      "MarkerFaceColor", "#FF8000")
plot3([x],[y],[0], ...
      'o', "Color",  "#FF8000", ...
      "MarkerFaceColor", "#FF8000")
plot3([x,x],[y,y],[0,z], ...
      '--', "Color",  "#FF8000")


xlabel("A1")
ylabel("A2")
zlabel("S1")

