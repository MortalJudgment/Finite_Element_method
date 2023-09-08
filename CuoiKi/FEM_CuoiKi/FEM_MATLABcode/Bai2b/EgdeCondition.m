function f = EgdeCondition(x,y,position)
switch position
    case 'Below'
        f = sin(2*pi*x) + x^2/4 + 1;
    case 'Left'
        f = 2*pi;
    case 'Above'
        f = sin(2*pi*x) + (x^2+1)/4 + 1;
    case 'Right'
        f = 2*pi + 1/2;
    otherwise
        disp('Not supported yet!');
end
