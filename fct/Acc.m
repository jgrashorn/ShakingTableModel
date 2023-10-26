function [ddx] = Acc(x,dx,k,kr,D,Dr,r0,rS0,rP0,m,I)

    % ACC Accelerations for the top plate,
    % x: state vector [3 positions, 3 rotations]
    % dx: velocity vector
    % k: 3x4 matrix containing x,y,z-stiffnesses of 4 springs
    % kr: 3x4 matrix of rotational stiffnesses
    % r0: 3x4 matrix containing vectors of spring positions (upper end)
    % rs0: 3x4 matrix containing vectors of spring positions (lower end)
    % rP0: 3x1 vector to COM of plate
    % m: plate mass
    % I: Moment of inertia of plate
    % D: translational damping
    % Dr: rotational damping

    rS = x(1:3)+Rotation(x)*r0; % actual position of the top of the springs
    rD = rS-r0; % displacement of the springs
    % rP = rS-rS0; % bottom of springs to top of springs
    rU = r0-rP0; % center of mass to top of springs
    rL = rP0-rS0+x(1:3); % center of mass to bottom of springs

    ddx_k = zeros(3,1);
    ddx_r = zeros(3,1);

    I_S = zeros(3,3,4);

    % Steiner
    for i=1:4
        A = [0,-rU(3,i), rU(2,i); rU(3,i), 0, -rU(1,i); -rU(2,i), rU(1,i), 0];
        I_S(:,:,i) = diag(I) + m * (A') * A;
    end

    ddx_k = zeros(3,4);
    ddx_r = zeros(3,4);

    for i=1:4
        F_r = Rotation(x)*-diag(k(:,i))*rD(:,i) - [0;0;m*9.81];
        F_d = Rotation(x)*-D*(dx(1:3)+cross(dx(4:6),rU(:,i)));
        ddx_k(:,i) = F_r./m + F_d./m; %- D*dx(1:3)/m;
        ddx_r(:,i) = (...
            - (cross(rL(:,i),F_r))./I... % Moments r \cross -kx
            - (cross(rL(:,i),F_d))./I... % Moments due to damping?
            - cross(dx(4:6),diag(I)*dx(4:6))./I... % Drall (?)
            - diag(kr(:,i))*x(4:6)./I... % restoring moments of springs
            - Dr*dx(4:6)./I); % damping
    end

    ddx = [sum(ddx_k,2);sum(ddx_r,2)];

end