function f = EgdeCondition(x,y,position)
switch position
    case 'Below'
        f = sin(2*pi*x) + x^2/4 + 1;
    case 'Left'
        f = cos(2*pi*y) + y^2/4 +2*pi;
    case 'Above'
        f = sin(2*pi*x) + (x^2+1)/4 + 1;
    case 'Right'
        f = cos(2*pi*y) + (y^2+1)/4 +2*pi + 1/2;
    otherwise
        disp('Not supported yet!');
end
