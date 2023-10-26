function [stateOut] = StateTransitionSimple(pf,state,dxdtF,t,meas)

    stateOut = zeros(size(state));
    nV = pf.NumStateVariables;

    % if length(unique(state(:,3))) < pf.NumParticles/2
    %     AIC = zeros(1,10);
    %     GMModels = cell(1,10);
    %     options = statset('MaxIter',2000);
    %     for k = 1:10
    %         try
    %             GMModels{k} = fitgmdist(state,k,'Options',options,'CovarianceType','diagonal');
    %             AIC(k)= GMModels{k}.AIC;
    %         catch
    %             AIC(k) = inf;
    %         end
    %     end
    % 
    %     [minAIC,numComponents] = min(AIC);
    %     if isinf(minAIC)
    %         warning("No GMM found!");
    %     else
    %         state = random(GMModels{numComponents},size(state,1));
    %     end
    % end

    parfor i=1:size(state,1)
        K = diag(state(i,5));
        D = diag(state(i,6));
        alpha = state(i,7);

        [~,stateOut_] = ode15s(dxdtF(K,D,alpha),t,state(i,1:2)');

        stateOut_ = stateOut_(end,:);
        stateOut(i,:) = [stateOut_,state(i,3:end).*(1+randn(size(state(i,3:end)))/100)];
    end

end