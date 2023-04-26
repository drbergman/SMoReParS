function [f_xy] = bilinearly_interpolate(X,Y,F,x,y)
% This function bilinearly interpolates a two-variable function F to find
% an estimate for F(x,y).

% find nearest elements of X,Y to the inputs x,y
index_x = 0;
index_y = 0;
found_x = false;
found_y = false;

for i= 1:length(X)-1
    if x <= X(i+1) && x >= X(i) && ~found_x
        index_x = i;
        found_x = true;
    end
    if y <= Y(i+1) && y >= Y(i) && ~found_y
        index_y = i;
        found_y = true;
    end
end

% bilinearly interpolate
dX = [X(index_x+1)- x, x - X(index_x)];
dY = [Y(index_y+1) - y; y - Y(index_y)];

a = 1/((X(index_x+1)-X(index_x))*(Y(index_y+1)-Y(index_y)));

sub_F = [F(index_x, index_y),  F(index_x,index_y+1);
         F(index_x+1,index_y), F(index_x+1,index_y+1)];

f_xy = a .* dX * sub_F * dY;
end

