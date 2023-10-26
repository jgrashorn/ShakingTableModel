function [lkl] = LogLikelihoodFull(theta,meas,dxdtFunc,t)

    bnd = [3500,.1,.6,1e-3;
           4500,1.5,1.4,1e-1];
    thetaPhys = bnd(1,:)+normcdf(theta).*(bnd(2,:)-bnd(1,:));

    pred = predFunc(thetaPhys,dxdtFunc,t,size(meas,1));

    lkl = zeros(size(pred,1),size(pred,2));

    for j=1:size(pred,2)
        for i=1:size(pred,3)
            lkl(:,j) = lkl(:,j) + log(normpdf(squeeze(pred(:,j,i)),meas(j,i),thetaPhys(:,4))+realmin);
        end
    end

    lkl = sum(lkl,2);

end

function[pred] = predFunc(theta,dxdtF,t,N)

    pred = zeros(size(theta,1),N);
    parfor i=1:size(theta,1)
        [~,pred_] = ode15s(dxdtF(theta(i,1),theta(i,2),theta(i,3)),t,zeros(2,1));
        pred(i,:) = pred_(:,1)';
    end

end
