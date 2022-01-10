% Finds the radius of a line segment starting at the origin
% and at an angle phi, and an ellipse with x-axis a, and y-axis b
%function r = intlineellipse(a,b,phi);
function r = intlineellipse(a,b,phi);

r = b*a/sqrt((a^2-b^2)*sin(phi)^2+b^2);
end
