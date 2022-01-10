%
%
%	function [g,k,s] = mrewind(k0,g0 [,gmax [,smax [,T]]])
%
%	Function tries to rewind the gradient optimally.
%
%	INPUT:
%		k0 = kx + i*ky of starting k-space point (cm^-1)
%		g0 = gx + i*gy of starting gradient (G/cm)
%		[gmax] = maximum gradient amplitude (default 3.9 G/cm)
%		[smax] = maximum slew rate (default 14500 G/cm/s)
%		T = time step (default 0.000004 s)
%
%	OUTPUT:
%		g = gradient waveform to return.
%		k = corresponding k-space trajectory.
%		s = corresponding slew rate.
%	


% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: rewind.m,v $
%	Revision 1.2  2002/08/16 23:30:22  bah
%	no message
%	
%	Revision 1.1  2002/03/28 01:07:28  bah
%	Added to cvs
%	
%
% ===========================================================


function [g,k,s] = mrewind(k0,g0,gmax,smax,T)


if (nargin < 3)
	gmax = 3.9;	% G/cm.
end;
if (nargin < 4)
	smax = 14500;	% G/cm/s.
end;
if (nargin < 5)
	T = .000004;	% seconds.
end;

gamma = 4258;		% Hz/G

gain = 20;	
plotstep = 5;	
kacc = .0001;		% Accuracy with which to rewind k.
gacc = 4/50000;		% Accuracy with which to rewind gradient.

rotang = angle(k0);	% First rotate trajectory.  This is to try to get
			% design to be insensitive to starting rotation.
k = abs(k0);
g = g0*exp(-i*rotang);

n = 1;
keepgoing = 1;

while (keepgoing==1)
    kret = k + gamma* 0.5 * g(n) * abs(g(n))/smax;	
				% k if we ramped g to 0 ASAP.
    kr(n) = kret;

    if (abs(kret) < kacc)	% Then keep ramping toward zero.
	ug = g(n)/abs(g(n));	% Unit vector along g.

	if (abs(g(n))>smax*T)
		g(n+1) = g(n) - ug * smax * T;
	else
		g(n+1) = 0;
	end;
	gtarg = 0;
    else
	gtarg = gain * (-kret);
	if (abs(gtarg)>gmax)
		gtarg = gmax * gtarg/abs(gtarg);
	end;
	dgvect = gtarg - g(n);
	if (abs(dgvect) > smax * T)
		dg = dgvect/abs(dgvect)*smax*T;
	else
		dg = dgvect;
	end;
        g(n+1) = g(n) + dg;
    end;
    karr(n) = k;

    if (1==0)	% ======= Animate Design ========= 
	
      if (rem(n,plotstep)==0)
	figure(1);
	plot(real(karr),imag(karr),'b-');
	hold on;
	plot(real(k),imag(k),'bx');
	plot(real([k kret]),imag([k kret]),'r-');
	plot(real([kret]),imag([kret]),'r+');
	title('k-space trajectory and projected return point');
	hold off;
	axis(1.1*abs(k)*[-1 1 -1 1]);
	axis square;
	grid on;

	figure(2);
	plot(real(g(1:length(g)-1)),imag(g(1:length(g)-1)),'b-');
	hold on;
	plot(real(g(n+1)),imag(g(n+1)),'bx');
	plot(real(gtarg),imag(gtarg),'r+');
	plotucirc(4);
	title('gradient trajectory and target');
	axis equal;
	hold off;
	grid on;
		
	pause;	
      end;
	
    end;


    % === ITERATE k... ====== %
	
%    k = k+gamma*(g(n)+g(n+1))/2*T;
    k = k+gamma*(g(n+1))*T;
    n = n+1;



 
    if (n > 3000)
	keepgoing = 0;
    end;
    if ((abs(k) < kacc) & (abs(g(n)) < gacc))
 	keepgoing = 0;
    end;
end;






%disp('   REWINDER GRADIENT DESIGN');
%disp('--------------------------------------------------');
%tt = sprintf('Final K-Space Amplitude = %f cm^(-1)',abs(k));
%disp(tt);
%tt = sprintf('Final Gradient Amplitude = %f G/cm',abs(g(n)));
%disp(tt);
%tt = sprintf('Rewinder Duration is %d us',round(1000000*T*length(g)));
%disp(tt);

k = [k0 karr]*exp(i*rotang);
g = g*exp(i*rotang);
s = ([g]-[g0 g(1:length(g)-1)])/T;






if (0==1)	% ======== Plot summary! ========= %
    
    figure(1);
    plot(real(k),imag(k));
    grid on;
    title('k-trajectory');
    axis equal;
    
    figure(2);
    plot(real(g),imag(g),'r-');
    hold on;
    plotucirc(gmax);
    title('Gradient Trajectory');
    
    t = [1:length(g)]*T;
    figure(3);
    subplot(2,2,1);
    plot(t,real(g),'b--',t,imag(g),'r--',t,abs(g),'b:');
    legend('Gx','Gy','Magnitude');
    title('Gradients vs Time');
    xlabel('Time(s)');
    ylabel('Gradient (G/cm)');
    grid on;
    
    subplot(2,2,2);
    plot(t,real(k),'b--',t,imag(k),'r--',t,abs(k),'b:');
    legend('Kx','Ky','Magnitude');
    title('K-vector vs Time');
    xlabel('Time(s)');
    ylabel('K (/cm)');
    grid on;
    
    subplot(2,2,3);
    plot(t,real(s),'b--',t,imag(s),'r--',t,abs(s),'b:');
    legend('Sx','Sy','Magnitude');
    title('Slew-vector vs Time');
    xlabel('Time(s)');
    ylabel('Slew Rate (G/cm/s)');
    grid on;
    
end;










