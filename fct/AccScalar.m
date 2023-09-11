function [ddx] = AccScalar(x,dx,k,kr,D,Dr,r0,rs0,rP0,m,I)

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

    rS = rs0 + x(1:3); % actual position of the plate's COM

    ddx_k = zeros(3,1);
    ddx_r = zeros(3,1);

    A = zeros(3,3,4);
    I_S = zeros(3,3,4);

    % Steiner
    for i=1:4
        A(:,:,i) = [0,-rS(3,i), rS(2,i)+x(2); rS(3,i)+x(3), 0, -rS(1,i); -rS(2,i), rS(1,i), 0];
        I_S(:,:,i) = diag(I) + m * (A') * A;
    end

    for i=1:4
        ddx_k = ddx_k - Rotation(-x)*diag(k(:,i))*(x(1:3))./m - D*dx(1:3)/m;
        ddx_r = ddx_r +(...
            - I_S(:,:,i)^-1 * (cross(rs0(:,i)+x(1:3),diag(-k(:,i))*Rotation(-x)*x(1:3)))... % Moments r \cross -kx
            - cross(dx(4:6),diag(I)*dx(4:6))./I... % Drall (?)
            - diag(kr(:,i))*x(4:6)./I... % restoring moments of springs
            - Dr*dx(4:6)./I); % damping
    end

    ddx = [ddx_k-[0;0;9.81];ddx_r];

end