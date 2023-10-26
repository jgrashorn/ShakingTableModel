function [stateOut] = StateTransition(pf,state,dxdtF,t)

    stateOut = zeros(size(state));
    nV = pf.NumStateVariables;

    % if size(unique(state(:,13))) < pf.NumParticles/10
    %     fprintf("Resampling...\n");
    %     state(:,13:end) = state(:,13:end).*(1+randn(size(state(:,13:end)))/10);
    % end

    if length(unique(state(:,3))) < pf.NumParticles/2
        AIC = zeros(1,10);
        GMModels = cell(1,10);
        options = statset('MaxIter',2000);
        for k = 1:10
            try
                GMModels{k} = fitgmdist(state,k,'Options',options,'CovarianceType','diagonal');
                AIC(k)= GMModels{k}.AIC;
            catch
                AIC(k) = inf;
            end
        end
        
        [minAIC,numComponents] = min(AIC);

        if isinf(minAIC)
            warning("Could not fit Gaussian model");
        else
            fprintf("Resampled with %d-dim Gaussian...\n",numComponents);
            state = random(GMModels{numComponents},size(state,1));
        end
    end

    parfor i=1:size(state,1)
        k = repmat(state(i,15:17)',1,4);
        k = state(i,21:23)'.*k;
        kr = state(i,21:23)'.*repmat(state(i,18:20)',1,4);

        D = state(i,20)*diag([1,1,1]); % lateral damping
        Dr = state(i,21)*diag([1,1,1]); % rotational damping

        [~,stateOut_] = ode15s(dxdtF(k,kr,D,Dr),t,state(i,1:12)');

        stateOut_ = stateOut_(end,:);
        stateOut(i,:) = [stateOut_,state(i,13:end).*(1+randn(size(state(i,13:end)))*.01)];
    end

    


end