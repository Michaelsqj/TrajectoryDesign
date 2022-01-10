%function [g] = gencone(RES,FOV,NINT,THETA,MAXLEN,TS,SMAX,GMAX,OVERSAMPLE,DCF);
% RES - Desired Resolution in mm
% FOV - Desired Field of View in cm
% NINT - Desired Number of Interleaves
% THETA - Range of polar angles (e.g. [0 pi/2] for 0 to 90 degrees)
% MAXLEN - Maximum length output allowed
% TS -  Sampling period in seconds
% SMAX - Slew rate (in G/cm/s)
% GMAX - Max Gradient Amplitude (in G)
% OVERSAMPLE - Amount to oversample 
% DCF - Density Compensation Factor

function [g] = gencone(RES,FOV,NINT,THETA,MAXLEN,TS,SMAX,GMAX,OVERSAMPLE,DCF);

if (nargin<5)
	MAXLEN = 10000;
end
if (nargin<6)
	TS = 0.000004;
end
if (nargin<7)
	SMAX = 15000;
end
if (nargin<8)
	GMAX = 3.98;
end
if (nargin<9)
	OVERSAMPLE = 4;
end
if (nargin<10)
	DCF = 1.0;
end
if (length(FOV)<2)
	FOV(2) = FOV(1);
end
if (length(RES)<2)
	RES(2) = RES(1);
end
THETA = min(abs(THETA),pi/2);
        RESrange = 1/(sqrt( (1/RES(1)*cos(min(abs(THETA))))^2+(1/RES(2)*sin(max(abs(THETA))))^2));
	THETArange = atan((1/RES(2)*sin(max(abs(THETA))))/(1/RES(1)*cos(min(abs(THETA)))));

FOVcirc = FOV(1);
maxtheta = max(abs(THETA));
mintheta = min(abs(THETA));
truemaxtheta = angle(cos(maxtheta)*RES(1)+i*sin(maxtheta)*RES(2));
truemintheta = angle(cos(mintheta)*RES(1)+i*sin(mintheta)*RES(2));

FOVrad = 1/max(intlineellipse(1/FOV(1),1/FOV(2),truemaxtheta),intlineellipse(1/FOV(1),1/FOV(2),truemintheta));
GTWISTADJUST = sqrt(1+sin(min(abs(THETA)))^2/max(0.0001,sin(max(abs(THETA))))^2*tan(min(abs(THETArange)))^2) / sqrt(1+tan(THETArange)^2); 

TS = TS/OVERSAMPLE;
[g,k,len] = wc(THETArange,[FOVcirc FOVrad],DCF,5/RESrange,NINT,GTWISTADJUST,MAXLEN*OVERSAMPLE,TS,SMAX*TS,GMAX);
TS = TS*OVERSAMPLE;
g = g(1:OVERSAMPLE:len,:);
% g(:,1) = smoothtraj(g(:,1),TS,SMAX);
% g(:,2) = smoothtraj(g(:,2),TS,SMAX);
% g(:,3) = smoothtraj(g(:,3),TS,SMAX);


end

function [gout] = smoothtraj(gin,Ts,smax);
des_sumg = cumsum(gin);

gout(1) = gin(1);
cursumg(1) = gout(1);

for (ai=(2:length(gin)))
	desg = des_sumg(ai)-cursumg(ai-1);
	dgin = desg-gout(ai-1);

	if (dgin>smax*Ts)
  	   gout(ai) = gout(ai-1)+smax*Ts;
        elseif (dgin<-smax*Ts)
           gout(ai) = gout(ai-1)-smax*Ts; 
        else
           gout(ai) = desg;
	end
	cursumg(ai) = cursumg(ai-1)+gout(ai); 
	if (abs(cursumg(ai) - des_sumg(ai))>0.00001)
		%disp('k-wise differences at');
		%ai
		%disp('of');
		%(cursumg(ai)-des_sumg(ai))*4258*Ts
	end
end

end