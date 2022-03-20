function idea_GradientResonance(g)

f1=1140;
f2=590;
bw1=220;
bw2=100;
oversample=1;
L=length(g)*oversample;
L=L-rem(L,2);
dt=(10e-6);

ft=fftshift(abs(fft(g,L)));
ft=ft(L/2+1:end);
f=linspace(0,1/dt/2,L/2);
figure,hold on

plot([f1-bw1/2 f1-bw1/2],[0,max(ft)],'.-r')
plot([f1+bw1/2 f1+bw1/2],[0,max(ft)],'.-r')
plot([f2-bw2/2 f2-bw2/2],[0,max(ft)],'.-r')
plot([f2+bw2/2 f2+bw2/2],[0,max(ft)],'.-r')

% axis([0 4000 0 1]);
plot(f,(ft),'k','LineWidth',1); xlim([0,2000]); ylim([0,max(ft)]);