function [tPhys] = ToPhysical(theta,bnd)
    
    tPhys = bnd(:,1) + normrnd(theta)*(bnd(2,:)-bnd(1,:));

end

