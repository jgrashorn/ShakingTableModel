function [R] = Rotation(x)
    
    Rx = @(x) [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
    Ry = @(x) [cos(x) 0 sin(x); 0 1 0; -sin(x) 0 cos(x)];
    Rz = @(x) [cos(x) -sin(x) 0; sin(x) cos(x) 0; 0 0 1];
    
    R = Rx(x(4))*Ry(x(5))*Rz(x(6));

end

