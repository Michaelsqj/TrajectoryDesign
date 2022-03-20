function [phi] = calc_phase(k)
    phi = atan2(k(:,2),k(:,1));
    phi(phi<0) = phi(phi<0)+2*pi;
    dphi = [0; diff(phi)]';
    dphi(dphi<0) = dphi(dphi<0)+ 2*pi;
    phi = cumsum(dphi);
end
