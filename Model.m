%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          ___________________________________________
%         /.                                          / 
%        /  A4 . . . . . . . . . . . . . . . . . . . / | l(3)
%       / .              z                          / A3
%      / .               |                         / /
%     / .               0| ___y                   / /
%    / .                /                        / /
%   / .              x /                    l(1)/ /
%  / .                                         / /
% / .                                         / /
% ____________________l(2)___________________/ /
% |                                         | /
% A1_______________________________________A2/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
addpath("fct/");

animate = 1; % use with caution
saveAnimation = 1;
fName = 'pics/Plate';

m = .210;
l = [75,75,15]; %distance to center of mass of corners
hSpring = 50; % height of springs, has no use currently

I = [1/12 * m * (l(2).^2 + l(3).^2); 1/12*m*(l(1).^2+l(3).^2); 1/12*m*(l(1).^2 + l(2).^2)];

stiffnessFactor = 25;
rStiffnessFactor = 2500;
kxNom = stiffnessFactor*10; kyNom = stiffnessFactor*15; kzNom = stiffnessFactor*15;
krxNom = rStiffnessFactor*10; kryNom = rStiffnessFactor*15; krzNom = rStiffnessFactor*12;

randomStiffness = 0;
rndFactor = .01;

% upper positions of springs
xA1_ = [l(1);-l(2);hSpring];
xA2_ = [l(1);l(2);hSpring];
xA3_ = [-l(1);l(2);hSpring];
xA4_ = [-l(1);-l(2);hSpring];

% lower positions of springs (in case its needed, idk)
xS1_ = xA1_ - [0;0;hSpring];
xS2_ = xA2_ - [0;0;hSpring];
xS3_ = xA3_ - [0;0;hSpring];
xS4_ = xA4_ - [0;0;hSpring];

switch randomStiffness
    case 0
        kA1x = kxNom; kA1y = kyNom; kA1z = kzNom;
        kA2x = kxNom; kA2y = kyNom; kA2z = kzNom;
        kA3x = kxNom; kA3y = kyNom; kA3z = kzNom;
        kA4x = kxNom; kA4y = kyNom; kA4z = kzNom;
    case 1
        kA1x = kxNom*(1+rndFactor*randn); kA1y = kyNom*(1+rndFactor*randn); kA1z = kzNom*(1+rndFactor*randn);
        kA2x = kxNom*(1+rndFactor*randn); kA2y = kyNom*(1+rndFactor*randn); kA2z = kzNom*(1+rndFactor*randn);
        kA3x = kxNom*(1+rndFactor*randn); kA3y = kyNom*(1+rndFactor*randn); kA3z = kzNom*(1+rndFactor*randn);
        kA4x = kxNom*(1+rndFactor*randn); kA4y = kyNom*(1+rndFactor*randn); kA4z = kzNom*(1+rndFactor*randn);
    otherwise
        error("Not implemented!");
end

k = [kA1x,kA2x,kA3x,kA4x;
     kA1y,kA2y,kA3y,kA4y;
     kA1z,kA2z,kA3z,kA4z];

kr = repmat([krxNom;kryNom;krzNom],1,4);

r0 = [xA1_,xA2_,xA3_,xA4_]; % positions of upper end of springs
rs0 = [xS1_,xS2_,xS3_,xS4_]; % positions of lower end of springs
rP0 = [0;0;hSpring+l(3)]; % position of COM of plate

% initial conditions
x0 = zeros(6,1); % position of the plate
x0(3) = -m*9.81/(sum(k(3,:))); % displacement due to gravity
dx0 = zeros(6,1); % velocity of the plate

% x0(1) = 10;
% x0(2) = 1;
% x0(3) = 1;
% x0(4) = 1*pi/180;
% x0(5) = 2*pi/180;
% x0(6) = 15*pi/180;

D = 0.1*eye(3); % lateral damping
Dr = 0.1*eye(3); % rotational damping

inputSwitch = 2; % 1: fake, 2: real

switch inputSwitch
    case 1
        %fake input
        T = 10;
        t = linspace(0,T,1000);
        % fx = @(t) 0*(5*sin(1.6*2*pi*t)+2*sin(4.3*2*pi*t)+sin(7*2*pi*t)).*(t<20);
        fx = @(t) 1.0 *sin(1*2*pi*t);
%         fx = @(t) 0 * t;
        fy = @(t) 0.0 * fx(t);
        f = @(t) [fx(t); fy(t); 0; 0; 0; 0];
    case 2
        %real input
        inp = load("input.mat","input");
        t = inp.input(:,1);
        fxReal = inp.input(:,2);
        fx = @(t_) interp1(t,fxReal,t_); % interpolating for ode-functions
        f = @(t_) [fx(t_); 0; 0; 0; 0; 0];
    otherwise
        error("Schrödingers input")
end

dxdt = @(t,x) [x(7:12);Acc((x(1:6)-f(t)),x(7:12),k,kr,r0,rs0,rP0,m,I,D,Dr)];

[t,x] = ode15s(dxdt,t,[x0;dx0]);

%% corner point coordinates
xA1 = zeros(3,length(t));
xA2 = xA1; xA3 = xA1; xA4 = xA1;
for i=1:length(t)
    xA1(:,i) = Rotation(x(i,:))*(x(i,1:3)' + r0(:,1));
    xA2(:,i) = Rotation(x(i,:))*(x(i,1:3)' + r0(:,2));
    xA3(:,i) = Rotation(x(i,:))*(x(i,1:3)' + r0(:,3));
    xA4(:,i) = Rotation(x(i,:))*(x(i,1:3)' + r0(:,4));
end

%% ANIMATION

framerate = 30;
imgId = 1;
oldT = 0;

if animate
    hFig = figure;
    for i=1:3
        plotLims{i} = [min([min(xA1(i,:)),min(xA2(i,:)),min(xA3(i,:)),min(xA4(i,:))]),...
                       max([max(xA1(i,:)),max(xA2(i,:)),max(xA3(i,:)),max(xA4(i,:))])];
    end
    if plotLims{3}(2)-plotLims{3}(1) < 5
        plotLimsDist = plotLims{3}(2)-plotLims{3}(1);
        plotLims{3}(1) = plotLims{3}(1)-(2.5-(plotLimsDist/2));
        plotLims{3}(2) = plotLims{3}(2)+(2.5-(plotLimsDist/2));
    end
    
    for i=1:length(t)
        xPlot = cell(3,1);
        for j=1:3
            xPlot{j} = [xA1(j,i),xA2(j,i),xA3(j,i),xA4(j,i)];
        end
        clf;
        
        xlim([plotLims{1}]);
        ylim([plotLims{2}]);
        zlim(plotLims{3});
    
        patch(xPlot{1},xPlot{2},xPlot{3},'b');
        view(3);
        title(num2str(t(i)));
        if saveAnimation
            if (t(i)-oldT)>1/framerate || imgId == 1
                exportgraphics(hFig, [fName  num2str(imgId) '.tif']);
                oldT = t(i);
                imgId = imgId + 1;
            end
        else
            pause(.01);
        end
    end
end

%% response plot
inp = zeros(6,length(t));
for i=1:length(t)
    inp(:,i) = f(t(i));
end
figure;
subplot(2,1,1); hold on;
for i=1:3
    plot(t,inp(i,:),'--');
    plot(t,x(:,i));
end
legend({"$f_x$","x","$f_y$","y","$f_z$","z"})
subplot(2,1,2); hold on;
for i=4:6
    plot(t,x(:,i));
end
legend({"$\theta_x$", "$\theta_y$","$\theta_z$"});

%% FFT plot
[ampX,fX] = fourier(x(:,1),t);
[ampY,fY] = fourier(x(:,2),t);
[ampZ,fZ] = fourier(x(:,3),t);
[ampIn,fIn] = fourier(fx(t),t);

figure;
plot(fX,ampX,fY,ampY,fZ,ampZ,fIn,ampIn);xlim([1,25]);
legend("X","Y","Z","Input");
title("FFT")
xlabel("Frequency [Hz]")
ylabel("Amplitude [mm]")
%% 3D line plot

figure;
subplot(2,4,1);
plot3(xA1(1,:),xA1(2,:),xA1(3,:));
subplot(2,4,2);
plot3(xA2(1,:),xA2(2,:),xA2(3,:));
subplot(2,4,5);
plot3(xA3(1,:),xA3(2,:),xA3(3,:));
subplot(2,4,6);
plot3(xA4(1,:),xA4(2,:),xA4(3,:));

subplot(2,4,[3,4,7,8])
hold on;
plot3(xA1(1,:),xA1(2,:),xA1(3,:));
plot3(xA2(1,:),xA2(2,:),xA2(3,:));
plot3(xA3(1,:),xA3(2,:),xA3(3,:));
plot3(xA4(1,:),xA4(2,:),xA4(3,:));
view(3)
xlabel("x");ylabel("y");zlabel("z");


%% Animation

function [ddx] = Acc(x,dx,k,kr,r0,rs0,rP0,m,I,D,Dr)

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

    rP = rP0 + x(1:3); % actual position of the plate's COM

    ddx_k = zeros(3,1);
    ddx_r = zeros(3,1);

    % Steiner
    A = [0,-rP(3), rP(2)+x(2); rP(3)+x(3), 0, -rP(1); -rP(2), rP(1), 0];
    I_S = diag(I) + m * (A') * A;

    for i=1:4
        ddx_k = ddx_k - Rotation(-x)*diag(k(:,i))*(x(1:3))./m - D*dx(1:3);
        ddx_r = ddx_r +(...
            - I_S^-1 * (cross(r0(:,i)+x(1:3),diag(k(:,i))*-Rotation(-x)*x(1:3)))... % Moments r \cross -kx
            - cross(dx(4:6),diag(I)*dx(4:6))./I... % Drall (?)
            - diag(kr(:,i))*x(4:6)./I... % restoring moments of springs
            - Dr*dx(4:6)); % damping
    end

    ddx = [ddx_k-[0;0;9.81];ddx_r];

end