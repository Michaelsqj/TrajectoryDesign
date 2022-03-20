function [Gs_new] = inPlaneRot(Gs, Gr, theta)
    % rotate Gs aroung Gr, because Gr defines the golden ratio direction
    % and is where the central axis lies
    Gs = Gs(:) ./ norm(Gs);
    Gr = Gr(:) ./ norm(Gr);
    Gp = cross(Gr, Gs);
    Gp = Gp ./ norm(Gp);
    
    % Gs_new is now on the plane defined by Gs, Gp
    % which are orthogonal to each other. and the angle
    % between Gs_new and Gs is theta, so now Gs_new can
    % be represented by the combination of Gs & Gp
    
    Gs_new = cos(theta).* Gs + sin(theta).*Gp;
    Gs_new = Gs_new ./ norm(Gs_new);
end