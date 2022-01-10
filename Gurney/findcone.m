%function [g,gr,nint] = findcone(RES,FOV,LEN,THETA,PRECISION,TS,WR,OS,SMAX,GMAX);
% g - Cone Gradient
% gr - Rewinder for Cone Gradient
% nint - Number of interleaves
% RES - Desired Resolution in mm
% FOV - Desired Field of View in cm
% LEN - Desired Length
% THETA - Range of polar angles (e.g. [0 pi/2] for 0 to 90 degrees)
% PRECISION - Precision of number of interleaves 
% TS -  Sampling period in seconds
% WR - Include Rewinder in length (1) or not (0)
% OS -  Oversampling Rate 
% SMAX - Slew rate (in G/cm/s)
% GMAX - Max Gradient Amplitude (in G)
% DCF - Density Compensation Factor
function [g,gr,nint] = findcone(RES,FOV,LEN,THETA,PRECISION,TS,WR,OS,SMAX,GMAX,DCF);

if (nargin<5)
	PRECISION = 0.1;
end
if (nargin<6)
        TS = 0.000004;
end
if (nargin<7)
	WR = 0;
end
if (nargin<8)
	OS = 1;
end
if (nargin<9)
        SMAX = 15000;
end
if (nargin<10)
        GMAX = 3.98;
end
if (nargin<11)
	DCF = 1.0;
end

NINTlo = 1;
NINThi = 100;
LENlo = 10000000;
withrew = WR;

curNINT = NINThi;
oldlength = 0;
MAXLEN = LEN+5;

g = gencone(RES,FOV,10000000,THETA,10000,TS,SMAX,GMAX,OS,DCF);
if (length(g)>LEN)
	disp('Sorry, but it is impossible to achieve that resolution in that length of time');
	g = [];
	gr = [];
	nint = 0;
else

g = gencone(RES,FOV,curNINT,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF);
	gc = (g(:,1)+i*g(:,2))*exp(-i*withrew*angle(g(end,1)+i*g(end,2)));
	gr = mrewind(sum(gc)*4258*TS,gc(end),GMAX,SMAX,TS);
	grz = mrewind(sum(g(:,3))*4258*TS,g(end,3),GMAX,SMAX,TS);
	gx = real(gc);
	gy = imag(gc);
	gz = g(:,3);

	gxr  = real(gr)'; 
	gyr =  imag(gr)';
	gzr = grz';

	if (length(grz)>length(gr))
		gxr(length(gzr)) = 0;
		gyr(length(gzr)) = 0;
        elseif (length(grz)<length(gr))
                gzr(length(gyr)) = 0;
	end

	g = [gx gy gz];
	gr = [gxr gyr gzr];
	if (withrew==1) leng =  length(g)+length(gr); else leng = length(g); end


LENhi = leng;
while ((leng>LEN) && (1))
	NINTlo = NINThi;
	NINThi = NINThi*2;
	curNINT = NINThi;
	LENhi = leng;
	oldlength = max(MAXLEN+1,leng);
        g = gencone(RES,FOV,curNINT,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF);
	gc = (g(:,1)+i*g(:,2))*exp(-withrew*i*angle(g(end,1)+i*g(end,2)));
	gr = mrewind(sum(gc)*4258*TS,gc(end),GMAX,SMAX,TS);
	grz = mrewind(sum(g(:,3))*4258*TS,g(end,3),GMAX,SMAX,TS);
	gx = real(gc);
	gy = imag(gc);
	gz = g(:,3);

	gxr  = real(gr)'; 
	gyr =  imag(gr)';
	gzr = grz';

	if (length(grz)>length(gr))
		gxr(length(gzr)) = 0;
		gyr(length(gzr)) = 0;
        elseif (length(grz)<length(gr))
                gzr(length(gyr)) = 0;
	end

	g = [gx gy gz];
	gr = [gxr gyr gzr];
	if (withrew==1) leng =  length(g)+length(gr); else leng = length(g); end

end


%if (leng==oldlength)
%	disp('Sorry, but it is impossible to achieve that resolution in that length of time');
%	g = [];
%	nint = 0;
%else
if (1)

%     while ((LENlo>LEN)
      while ((((NINThi-NINTlo)>PRECISION) || (leng>LEN)) && (LENhi~=LEN) )
	curNINT = (NINThi-NINTlo)/2+NINTlo;
        g = gencone(RES,FOV,curNINT,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF);
	gc = (g(:,1)+i*g(:,2))*exp(-withrew*i*angle(g(end,1)+i*g(end,2)));
	gr = mrewind(sum(gc)*4258*TS,gc(end),GMAX,SMAX,TS);
	grz = mrewind(sum(g(:,3))*4258*TS,g(end,3),GMAX,SMAX,TS);
	gr = gr(2:end);
	grz = grz(2:end);
	gx = real(gc);
	gy = imag(gc);
	gz = g(:,3);

	gxr  = real(gr)'; 
	gyr =  imag(gr)';
	gzr = grz';

	if (length(grz)>length(gr))
		gxr(length(gzr)) = 0;
		gyr(length(gzr)) = 0;
        elseif (length(grz)<length(gr))
                gzr(length(gyr)) = 0;
	end

	g = [gx gy gz];
	gr = [gxr gyr gzr];
	if (withrew==1) leng =  length(g)+length(gr); else leng = length(g); end

	if (leng>LEN)
             NINTlo = curNINT;
	     LENlo = leng;
	else
             NINThi = curNINT;
	     LENhi = leng;
	end
	%[NINTlo NINThi LENlo LENhi leng]
	%plot(g);
	drawnow;
     end
end
nint = NINThi;
if (withrew ==1)
if (length(g)+length(gr))>LEN
	disp('Major ERROR');
	length(g)+length(gr)
end
else 
if (length(g))>LEN
	disp('Major ERROR');
	length(g)
end
end
if (withrew==1)
	if ((length(g)+length(gr))<LEN)
		gr(LEN-length(g),:)=0;
	end
end

end