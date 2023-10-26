function [out] = LikelihoodSimple(pf,state,msm,f,v,ind)

    out = ones(size(state,1),1);
    stateInd = state(:,3:end)<0;

    meas = msm';
    pred = state(:,ind);
    
    for i=1:size(msm,2)
        if ind(i)<2
            out = out.*normpdf(pred(:,i),meas(i,:),state(:,3)).*(sum(stateInd,2)==0);
        else
            out = out.*normpdf((pred(:,i)-meas(i,:)')/meas(i,:),state(:,4)).*(sum(stateInd,2)==0);
        end
    end

    out(isnan(out)) = 0;

end

