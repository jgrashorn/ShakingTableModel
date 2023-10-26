clear; close all;
addpath("fct/");

%% constants (hopefully)
m = .240;
l = [75,75,15]*1e-3; %distance to center of mass of corners
hSpring = 100e-3; % height of springs

I = [1/12 * m * (l(2).^2 + l(3).^2); 1/12*m*(l(1).^2+l(3).^2); 1/12*m*(l(1).^2 + l(2).^2)];

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

r0 = [xA1_,xA2_,xA3_,xA4_]; % positions of upper end of springs
rs0 = [xS1_,xS2_,xS3_,xS4_]; % positions of lower end of springs
rP0 = [0;0;hSpring+l(3)]; % position of COM of plate

% initial conditions
x0 = zeros(6,1); % position of the plate
dx0 = zeros(6,1); % velocity of the plate

%% data treatment

%% full data from before 29.09.2023
% scaling of the data was unsure, therefore manual scaling of 1.4 was used.

nDofs = 6;

X = readmatrix("../SINDy/test3/3DOF/X.csv");
dX = readmatrix("../SINDy/test3/3DOF/dX.csv");
Y = readmatrix("../SINDy/test3/3DOF/Y.csv");
dY = readmatrix("../SINDy/test3/3DOF/dY.csv");
tSim = readmatrix("../SINDy/test3/3DOF/t.csv");

X = X(:,1:nDofs);
Y = Y(:,1:nDofs);
dX = dX(:,1:nDofs);
dY = dY(:,1:nDofs);

xMeas = [X,dX];
yMeas = [Y,dY];

tStart = 1.76;
tStartInd = find(tSim>tStart,1);

for i=1:size(X,2)
    fxFunc{i} = @(t) interp1(tSim,Y(:,i),t);
    vxFunc{i} = @(t) interp1(tSim,dY(:,i),t);
    fxFuncNoise{i} = @(t) interp1(tSim,Y(:,i),t)*(1+randn*.01);
    vxFuncNoise{i} = @(t) interp1(tSim,dY(:,i),t)*(1+randn*.01);
end

if size(X,2)==3
    f = @(t) [fxFunc{1}(t);fxFunc{2}(t);fxFunc{3}(t)];
    v = @(t) [vxFunc{1}(t);vxFunc{2}(t);vxFunc{3}(t)];
elseif size(X,2)==6
    f = @(t) [fxFunc{1}(t);fxFunc{2}(t);fxFunc{3}(t);fxFunc{4}(t);fxFunc{5}(t);fxFunc{6}(t)];
    v = @(t) [vxFunc{1}(t);vxFunc{2}(t);vxFunc{3}(t);vxFunc{4}(t);vxFunc{5}(t);vxFunc{6}(t)];
else
    error("X dimensions are weird.");
end

measInd = [1,2,3,4,5,6]; % which indices of x are measured

%% function setup

simple = 1;

if ~simple
    dxdtPF = @(k,kr,D,Dr) @(t,x) [x(7:12);Acc((x(1:6)-f(t)),(x(7:12)-v(t)),k,kr,D,Dr,r0,rP0,m,I)];
else
    f = @(t) [fxFunc(t)];
    v = @(t) 0;
    measInd = [1,2];
    % xMeas = [xMeas(1,:);xMeas(2,:);xMeas(3,:)];
    dxdtPF = @(K,D,alpha) @(t,x) [x(2);-K*(alpha*x(1)-f(t))-D*(alpha*x(2)-vxFunc(t))];
end

%% complex particle filter initialization

if ~simple
    pf = stateEstimatorPF;
    pf.StateEstimationMethod = 'mean';
    pf.ResamplingMethod = 'systematic';
    
    pf.StateTransitionFcn = @StateTransition;
    pf.MeasurementLikelihoodFcn = @Likelihood;
    
    nParticles = 1000;
    
    initState = [x0; % [x0] 6x1, 1:6
                 dx0; % [dx0] 6x1, 7:12
                 1; % sigma pos, 13
                 10; % sigma vel, 14
                 1000*ones(3,1); % kd, 15:17
                 10*ones(3,1); % kr, 18:20
                 .0001; % D, 20
                 10; % Dr, 21
                 1; % mpl X, 22
                 1; % mpl Y, 23
                 1; % mpl Z, 24
                 ]';
    
    pfCOV = [0.01*ones(12,1); %x0, dx0, 1:12
             .1; %sigma pos, 13
             1; % sigma vel, 14
             200*ones(3,1); % kd,15:17
             1*ones(3,1); % kr, 18:20
             .00001; % D, 20
             1; % Dr, 21
             .1; % mpl X, 22
             .1; % mpl Y, 23
             .1 % mpl Z, 24
             ].^2;
    pfCOV = diag(pfCOV);
    
    initialize(pf, nParticles, initState, pfCOV);
    
    statePred = zeros(length(tSim),pf.NumStateVariables);
    statePred(1,:) = initState';
    stateCorrected = zeros(length(tSim),pf.NumStateVariables);
    
    figure;
    clrs = {'b','r','k','g','y','m'};
    
    for i=tStartInd:length(tSim)
            
        [statePred(i,:),covPred] = predict(pf,dxdtPF,[tSim(i-1),tSim(i)]);
    
        [stateCorrected(i,:), covCorrected] = correct(pf, xMeas(i,measInd)',f(tSim(i)),v(tSim(i)),measInd);
    
        if mod(i,5)==0
            fprintf("State %d of %d\n",i,length(tSim));
            clf;
            posInd = 1;
            velInd = 1;
            for j=measInd
                clear plotState
                for ii=tStartInd:i
                    if j<7
                        f_ = f(tSim(ii));
                        plotState(j,ii-(tStartInd-1)) = stateCorrected(ii,j);
                    else
                        v_ = v(tSim(ii));
                        plotState(j,ii-(tStartInd-1)) = stateCorrected(ii,j);
                    end
                end
                if j<7
                    subplot(2,2,1);hold on;
                    plot(tSim(tStartInd:i),plotState(j,:),[clrs{posInd} '-']);
                    plot(tSim(tStartInd:i),xMeas(tStartInd:i,j),[clrs{posInd} '--']);
                    posInd = posInd+1;
                elseif j<13
                    subplot(2,2,2);hold on;
                    plot(tSim(tStartInd:i),plotState(j,:),[clrs{velInd} '-']);
                    plot(tSim(tStartInd:i),xMeas(tStartInd:i,j),[clrs{velInd} '--']);
                    velInd = velInd+1;
                end
                % plot(t(1:i),measurementNoise(1:i,j));
            end
            legend("x","",'Location','southwest');
            ylim([min(min(xMeas(tStartInd:i,measInd))),max(max(xMeas(tStartInd:i,measInd)))]);    
            subplot(2,2,[3,4]);
            hold on;
            for j=1:pf.NumStateVariables-12
                plot(tSim(tStartInd:i),statePred(tStartInd:i,12+j)./initState(12+j))
            end
            drawnow;
        end
    
    end
%% 
else
    pf = stateEstimatorPF;
    pf.StateEstimationMethod = 'mean';
    pf.ResamplingMethod = 'systematic';
    
    pf.StateTransitionFcn = @StateTransitionSimple;
    pf.MeasurementLikelihoodFcn = @LikelihoodSimple;
    
    nParticles = 2000;

    stepParticles = zeros(nParticles,4,length(tSim));
    stepWeights = zeros(nParticles,length(tSim));
    
    initState = [zeros(2,1);.1;2;4000;.9;1]'; % x0,dx0,sigma_pos,sigma_vel,k,d,alpha
    
    pfCOV = [0.01*ones(2,1);.01;.1;600;.2;.2].^2;
    pfCOV = diag(pfCOV);
    
    initialize(pf, nParticles, initState, pfCOV);

    pfPolicy = resamplingPolicyPF;
    pfPolicy.TriggerMethod = 'interval';

    pf.ResamplingPolicy = pfPolicy;
    
    statePred = zeros(length(tSim),pf.NumStateVariables);
    statePred(1,:) = initState';
    stateCorrected = zeros(length(tSim),pf.NumStateVariables);
    
    figure;
    clrs = {'b','r','k'};

    measEvery = 3;
    
    for i=tStartInd:length(tSim)
        
        [statePred(i,:),covPred] = predict(pf,dxdtPF,[tSim(i-1),tSim(i)]);
        
        if i==tStart || mod(i,measEvery)==0
            [stateCorrected(i,:), covCorrected] = correct(pf, xMeas(measInd,i)',f(tSim(i)),v(tSim(i)),measInd);
        end

        stepParticles(:,1:2,i) = pf.Particles(:,1:2);
        stepParticles(:,3:5,i) = pf.Particles(:,5:7);
        stepWeights(:,i) = pf.Weights;
    
        if mod(i,5)==0
            fprintf("State %d of %d\n",i,length(tSim));
            clf;
            clrInd = 1;
            for j=measInd
                clear plotState
                for ii=tStartInd:i
                    if j<2
                        f_ = f(tSim(ii));
                        plotState(j,ii-(tStartInd-1)) = statePred(ii,j);
                    elseif j<3
                        v_ = v(tSim(ii));
                        plotState(j,ii-(tStartInd-1)) = statePred(ii,j);
                    end
                end
                if j<2
                    subplot(2,2,1);hold on;
                    plot(tSim(tStartInd:i),plotState(j,:),[clrs{clrInd} '-']);
                    plot(tSim(tStartInd:i),xMeas(j,tStartInd:i),[clrs{clrInd} '--']);
                elseif j<3
                    subplot(2,2,2);hold on;
                    plot(tSim(tStartInd:i),plotState(j,:),[clrs{clrInd} '-']);
                    plot(tSim(tStartInd:i),xMeas(j,tStartInd:i),[clrs{clrInd} '--']);
                end
                clrInd = clrInd+1;
                % plot(t(1:i),measurementNoise(1:i,j));
            end
            legend("x","",'Location','southwest');
            ylim([min(min(xMeas(measInd,tStartInd:i))),max(max(xMeas(measInd,tStartInd:i)))]);    
            subplot(2,2,[3,4]);
            hold on;
            for j=1:pf.NumStateVariables-2
                plot(tSim(tStartInd:i),statePred(tStartInd:i,2+j)./initState(2+j))
            end
            drawnow;
        end
    
    end
end

%%
for i=1:length(tSim)
    for j=1:5
        meanParticles(i,j) = mean(stepParticles(:,j,i));
        upperBndParticles(i,j) = meanParticles(i,j)+1.96*std(stepParticles(:,j,i));
        lowerBndParticles(i,j) = meanParticles(i,j)-1.96*std(stepParticles(:,j,i));
        minBndParticles(i,j) = min(stepParticles(:,j,i));
        maxBndParticles(i,j) = max(stepParticles(:,j,i));
    end
end

figure;
titles = {"$x$","$\dot{x}$","$k$","$d$","$\alpha$"};
for i=1:5
    subplot(2,3,i);
    plot(tSim,meanParticles(:,i),'r');hold on;
    fill([tSim,fliplr(tSim)],[upperBndParticles(:,i)',fliplr(lowerBndParticles(:,i)')],'b','FaceAlpha',.8,'EdgeAlpha',0);
    fill([tSim,fliplr(tSim)],[maxBndParticles(:,i)',fliplr(minBndParticles(:,i)')],'b','FaceAlpha',.3,'EdgeAlpha',0);
    plot(tSim,meanParticles(:,i),'r');
    if i<3
        plot(tSim,xMeas(i,:),'k--','LineWidth',2);
    end
    xlabel("Time [s]");
    ylabel(titles{i});
    xlim([tSim(tStartInd),tSim(tStartInd)+2])
    % set(gca,'Color','k')
    hl = legend({"Mean","95\% interval"});
    % set(hl, 'TextColor','w', 'Color','none', 'EdgeColor','none')
end
set(findobj(gcf,'type','axes'),'FontSize',25);

%% FEMU

acc = diff(diff(xMeas(1,:)))./(tSim(2)-tSim(1));
vel = vx;

for i=1:length(tSim)-2
    minK(i,:) = fminunc( @(in) norm(unbF(in,acc,vel(1:end-2),xMeas(1,1:end-2)-fx(1:end-2),i)),[4000,1]',optimoptions(@fminunc,'Display','none'));
end

function[out] = unbF(in,acc,vel,resp,i)
    
    dK = in(1);
    D = in(2);
    dim = size(acc,1);
    out = (acc(i)' + dK.*eye(dim)*(resp(i)') + D.*eye(dim)*vel(i)').^2;

end
