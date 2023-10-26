clear;

X = readmatrix("../SINDy/test3/3DOF/X.csv");
Y = readmatrix("../SINDy/test3/3DOF/Y.csv");
tSim = readmatrix("../SINDy/test3/3DOF/t.csv");

noiseIn = .01;
noiseOut = .05;

Y = Y(:,1).*(1+randn(size(Y(:,1)))*noiseIn);
dY(1) = 0;
dY(2:1700) = diff(Y)./(tSim(2)-tSim(1));

u1Func = @(t) interp1(tSim,Y,t);
u2Func = @(t) interp1(tSim,dY,t);

model = @(k,d) @(t,x) [x(2);-k*(x(1)-u1Func(t))-d*(x(2)-u2Func(t))];

k = 4000;
d = 1;

[~,real] = ode15s(model(4000,1),tSim,zeros(2,1));

plot(tSim,real(:,1))

realNoise = real.*(1+randn(size(real))*noiseOut);

plot(tSim,realNoise(:,1));

tmcmcSamples = supermethods.tmcmc_jan(@(t) LogLikelihoodFull(t,realNoise,model,tSim),...
                    @(x) mvnpdf(x),@(N) mvnrnd(ones(2,1),eye(2),N),1000);


%%
smpPhys = Transf(tmcmcSamples);
smpPhys = mean(smpPhys);

[~,updated] = ode15s(model(smpPhys(1),smpPhys(2)),tSim,zeros(2,1));

diffU = updated(:,1)-realNoise(:,1);
diffReal = real(:,1)-realNoise(:,1);

figure;
plot(tSim,realNoise(:,1),tSim,updated(:,1));
figure;
plot(tSim,diffU);

mean(diffU)
std(diffU)

mean(diffReal)
std(diffReal)


%% functions

function [phys] = Transf(t)

    bnd = [500,.1,;
           5000,1.5];

    phys = bnd(1,:)+normcdf(t).*(bnd(2,:)-bnd(1,:));

end
function [lkl] = LogLikelihoodFull(theta,meas,dxdtFunc,t)

   
    thetaPhys = Transf(theta);

    pred = predFunc(thetaPhys,dxdtFunc,t,size(meas,1));

    lkl = zeros(size(pred,1),size(pred,2));

    for j=1:size(pred,2)
        for i=1:size(pred,3)
            lkl(:,j) = lkl(:,j) + log(normpdf(squeeze(pred(:,j,i)),meas(j,i),1e-2)+realmin);
        end
    end

    lkl = sum(lkl,2);

end

function[pred] = predFunc(theta,dxdtF,t,N)

    pred = zeros(size(theta,1),N);
    parfor i=1:size(theta,1)
        [~,pred_] = ode15s(dxdtF(theta(i,1),theta(i,2)),t,zeros(2,1));
        pred(i,:) = pred_(:,1)';
    end

end