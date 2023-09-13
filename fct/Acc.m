function [ddx] = Acc(x,dx,k,kr,D,Dr,rS0,rP0,m,I)

    % Accelerations for the top plate,
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

    rS = Rotation(x)*(x(1:3)+rS0); % actual position of the top of the springs
    rD = rS-rS0;
    rP = rP0-rS;

    ddx_k = zeros(3,1);
    ddx_r = zeros(3,1);

    I_S = zeros(3,3,4);

    % Steiner
    for i=1:4
        A = [0,-rP(3,i), rP(2,i); rP(3,i), 0, -rP(1,i); -rP(2,i), rP(1,i), 0];
        I_S(:,:,i) = diag(I) + m * (A') * A;
    end

    for i=1:4
        F_r = diag(-k(:,i))*rD(:,i);
        ddx_k = ddx_k + F_r./m - D*dx(1:3)/m;
        ddx_r = ddx_r +(...
            - I_S(:,:,i)^-1 * (cross(rS(:,i),F_r))... % Moments r \cross -kx
            - cross(dx(4:6),diag(I)*dx(4:6))./I... % Drall (?)
            - diag(kr(:,i))*x(4:6)./I... % restoring moments of springs
            - Dr*dx(4:6)./I); % damping
    end

    ddx = [ddx_k-[0;0;0];ddx_r];

end