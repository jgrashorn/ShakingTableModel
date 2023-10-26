function [out] = Likelihood(pf,state,msm,f,v,ind)

    out = zeros(size(state,1),1);
    stateInd = state(:,13:end)<0;

    meas = msm;
    pred = state(:,ind);

    for i=1:size(meas,1)
        if ind(i)<7
            out = out + log(normpdf((pred(:,i)-meas(i,:))./meas(i,:),0,state(:,13)).*(sum(stateInd,2)==0));
        else
            out = out + log(normpdf((pred(:,i)-meas(i,:))./meas(i,:),0,state(:,14)).*(sum(stateInd,2)==0));
        end
    end
    out = exp(out);
    out(isnan(out)) = 0;

end

