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

animate = 0; % use with caution
saveAnimation = 0;
fName = 'pics/Plate';

if animate && saveAnimation && ~strcmp(input("You sure you want to save the animation? y/n ",'s'),'y')
    error("Canceled")
end

useParticleFilter = 0;
teapot = 0;

m = .240;
l = [75,75,15]*1e-3; %distance to center of mass of corners
hSpring = 100e-3; % height of springs

I = [1/12 * m * (l(2).^2 + l(3).^2); 1/12*m*(l(1).^2+l(3).^2); 1/12*m*(l(1).^2 + l(2).^2)];

% ratios of spring stiffnesses
xRatio = 1;
yRatio = 2;
zRatio = .8;

% rotational
xRRatio = 1;
yRRatio = 1;
zRRatio = 1;

stiffnessFactor = 20;
rStiffnessFactor = .5;
kxNom = stiffnessFactor*xRatio; kyNom = stiffnessFactor*yRatio; kzNom = stiffnessFactor*zRatio;
krxNom = rStiffnessFactor*xRRatio; kryNom = rStiffnessFactor*yRRatio; krzNom = rStiffnessFactor*zRRatio;

randomStiffness = 0;
rndFactor = .02;

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
% x0(3) = -m*9.81/(sum(k(3,:))); % displacement due to gravity
dx0 = zeros(6,1); % velocity of the plate

% x0(1) = -.01;
% x0(2) = -.01;
% x0(3) = -.01;
% x0(4) = 10*pi/180;
% x0(5) = -2*pi/180;
% x0(6) = 2*pi/180;

D = .0001*diag([1,1,1]); % lateral damping
Dr = .00001*diag([1,1,1]); % rotational damping

inputSwitch = 1; % 1: fake, 2: real

switch inputSwitch
    case 1
        %fake input
        T = 20;
        t = linspace(0,T,1000);
        fx = @(t) 1* (.005*sin(1.6*2*pi*t)+.002*sin(4.3*2*pi*t)+.001*sin(7*2*pi*t)).*(t<20);
        % fx = @(t) 8e-3 * sin(12*2*pi*t);
%         fx = @(t) 0 * t;
        fy = @(t) 0.0 * fx(t);
        f = @(t) [fx(t); fy(t); 0; 0; 0; 0];
    case 2
        %real input
        inp = load("input.mat","input");
        t = inp.input(:,1);
        dt = t(2)-t(1);
        tAdd = (t(end)+dt:dt:floor(t(end)*1.2));
        t = [t;tAdd'];
        fxReal = inp.input(:,2)*1e-3;
        fxReal = [fxReal;zeros(length(tAdd),1)];
        fx = @(t_) interp1(t,fxReal,t_); % interpolating for ode-functions
        f = @(t_) [fx(t_); 0; 0; 0; 0; 0];
    otherwise
        error("Schrödingers input")
end

dxdt = @(t,x) [x(7:12);Acc((x(1:6)-f(t)),x(7:12),k,kr,D,Dr,r0,rP0,m,I)];

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

framerate = 100;
imgId = 1;
oldT = 0;
minZ = .1;

if animate

    if teapot
        [verts_, faces, cindex] = teapotGeometry;
        verts = zeros([size(verts_),length(t)]);
        
        for i=1:length(t)
            verts(:,:,i) = (Rotation(x(i,:))*(verts_'+x(i,1:3)'*1e3))';
        end
    
        xV = verts(:,1,:);
        xL = [min(xV(:)),max(xV(:))];
        yV = verts(:,2,:);
        yL = [min(yV(:)),max(yV(:))];
        zV = verts(:,3,:);
        zL = [min(zV(:)),max(zV(:))];
    end

    hFig = figure;
    for i=1:3
        plotLims{i} = [min([min(xA1(i,:)),min(xA2(i,:)),min(xA3(i,:)),min(xA4(i,:))]),...
                       max([max(xA1(i,:)),max(xA2(i,:)),max(xA3(i,:)),max(xA4(i,:))])];
    end
    if plotLims{3}(2)-plotLims{3}(1) < minZ
        plotLimsDist = plotLims{3}(2)-plotLims{3}(1);
        plotLims{3}(1) = plotLims{3}(1)-(minZ/2-(plotLimsDist/2));
        plotLims{3}(2) = plotLims{3}(2)+(minZ/2-(plotLimsDist/2));
    end
    
    for i=1:length(t)
        if teapot
            clf;
            patch('Faces',faces,'Vertices',squeeze(verts(:,:,i)),'FaceVertexCData',cindex,'FaceColor','interp');
            view(3);
            ylim(xL);
            xlim(yL);
            zlim(zL);
        else
            xPlot = cell(3,1);
            pPlot = cell(3,1);
            fPlot = f(t(i));
            for j=1:3
                xPlot{j} = [xA1(j,i),xA2(j,i),xA3(j,i),xA4(j,i)];
                pPlot_ = zeros(1,4);
                if j==3
                    pPlot{j} = repmat(plotLims{3}(1),1,4);
                else
                    for jj=1:4
                        pPlot_(jj) = rs0(j,jj)+fPlot(j);
                    end
                    pPlot{j} = pPlot_;
                end
            end
            clf;
            
            xlim([plotLims{1}]);
            ylim([plotLims{2}]);
            zlim(plotLims{3});
            hold on;
            patch(xPlot{1},xPlot{2},xPlot{3},'m');
            patch(pPlot{1},pPlot{2},pPlot{3},'g');
            hold off;
            view(3);
            title(num2str(t(i)));
        end

        if saveAnimation
            if (t(i)-oldT)>1/framerate || imgId == 1
                exportgraphics(hFig, [fName  num2str(imgId) '.tif']);
                oldT = t(i);
                imgId = imgId + 1;
            end
        else
            pause(.005);
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

%% response of corner A1 in x-direction

try
    tData = load("testingData.mat");
    tMeas = tData.tMeas(tData.tMeas>1.41)-1.41;
    dataMeas = -tData.xMeas(tData.tMeas>1.41)*1e-3;
    figure; hold on;
    plot(t,xA1(1,:)-inp(1,:)-r0(1,1));
    plot(tMeas,dataMeas);
    legend("Simulated","Measured");

    [ampA1,fA1] = fourier(xA1(1,:)-r0(1,1),t);
    [ampMeas,fMeas] = fourier(dataMeas,tMeas);

    figure; hold on;
    plot(fA1,ampA1);
    plot(fMeas,ampMeas);
    xlim([0,25])
    legend("Simulated","Measured");
    xlabel("Frequency [Hz]");ylabel("Amplitude");
catch
    warning("Test data not found, continuing...");
end

%% Particle filter

if useParticleFilter
    close all
    
    mult0 = [.5,.5,.5]';
    dxdtPF = @(mpl) @(t,x) [x(7:12);Acc((x(1:6)-f(t)),x(7:12),mpl.*k,kr,D,Dr,r0,rP0,m,I)];

    [~,xM] = ode15s(dxdtPF(mult0),t,[x0;dx0]);
    
    pf = stateEstimatorPF;
    pf.StateEstimationMethod = 'mean';
    pf.ResamplingMethod = 'systematic';
    
    pf.StateTransitionFcn = @StateTransition;
    pf.MeasurementLikelihoodFcn = @Likelihood;
    
    nParticles = 2500;

    initState = [x0;dx0;.3;1;1;1]';
    
    initialize(pf, nParticles, initState, diag([0.0*ones(12,1);.02;.2;.2;.2]).*eye(16));

    statePred = zeros(length(t),pf.NumStateVariables);
    statePred(1,:) = initState';
    stateCorrected = zeros(length(t),pf.NumStateVariables);

    inp = zeros(6,length(t));
    for i=1:length(t)
        inp(:,i) = f(t(i));
    end

    measurement = xM(:,1:3)';
    % measurement = measurement;
    measurementNoise = measurement.*(1+randn(size(measurement))/5);
    % updateData = inpFile.xMeas(inpFile.tMeas>t_Start);
    updateData = measurementNoise-inp(1:3,:);
    
    toPlotIdx = randperm(nParticles,50);
    figure;
    clrs = {'b','r','k'};
    for i=2:length(t)
        
        if i==1
            [statePred(i,:),covPred] = predict(pf,dxdtPF,[0,t(i)]);
        else
            [statePred(i,:),covPred] = predict(pf,dxdtPF,[t(i-1),t(i)]);
        end
        [stateCorrected(i,:), covCorrected] = correct(pf, updateData(:,i)',f(t(i)));
    
        if mod(i,5)==0
            fprintf("State %d of %d\n",i,length(t));
            clf;
            subplot(2,1,1);
            hold on;
            for j=1:size(updateData,1)
                clear plotState
                for ii=1:i
                    plotState(:,ii) = Rotation(stateCorrected(ii,1:6))*stateCorrected(ii,1:3)'-inp(j,ii);
                end
                plot(t(1:i),plotState(j,:),[clrs{j} '-']);
                plot(t(1:i),updateData(j,1:i),[clrs{j} '--']);
                % plot(t(1:i),measurementNoise(1:i,j));
            end
            legend("x","","y","","z","",'Location','southwest');
            ylim([min(min(updateData(:,1:i))),max(max(updateData(:,1:i)))]);    
            subplot(2,1,2);
            hold on;
            for j=1:pf.NumStateVariables-12
                plot(t(1:i),statePred(1:i,12+j))
            end
            drawnow;
        end
    
    end
end

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

%% functions

function [stateOut] = StateTransition(pf,state,dxdtF,t)

    stateOut = zeros(size(state));
    nV = pf.NumStateVariables;
    parfor i=1:size(state,1)
        kMpl = state(i,14:16)';
    
        dxdt = dxdtF(kMpl);
    
        [~,stateOut_] = ode15s(dxdt,t,state(i,1:12)');
        stateOut_ = stateOut_(end,:);
        stateOut(i,:) = [stateOut_,state(i,13:nV)];
    end

end

function [out] = Likelihood(pf,state,msm,f)

    % if any(state(:,13:pf.NumStateVariables)<1e-6)
    %     out = zeros(size(state,1),1);
    %     return
    % end

    out = ones(size(state,1),1);
    stateInd = state(:,13:16)<1e-6;

    % meas = msm'-f(1:size(msm,2));
    % pred = state(:,1:3)-f(1:3)';

    meas = msm';
    pred = state(:,1:3)-f(1:3)';
    
    for i=1:size(msm,1)
        out = out.*normpdf(pred(:,i),meas(i,:)',state(:,13)).*(sum(stateInd,2)==0);
    end

    out(isnan(out)) = 0;

end