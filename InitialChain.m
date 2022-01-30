function [R0,N,r0] = InitialChain(Ns,b)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N=sum(Ns);
R0=RandomWalk(N, 3)*b;
w=wgn(N,3,0)*b/5;
R0=R0+w;
maxR=max(R0);minR=min(R0);offsetR=(maxR+minR)/2;
R0=R0-offsetR;
r2=max(sum(R0.^2,2));
%figure(),plot3(R0(:,1),R0(:,2),R0(:,3))
PlotRouseChain(R0,Ns,0);
r0=sqrt(r2);
x=(-1:0.01:1);
y=sqrt(1-x.^2);
x=[x,x(end:-1:1)];y=[y,-y];
x=x*r0;y=y*r0;
plot(x,y,'-','LineWidth',3,'MarkerSize',12,'Color','#77AC30')
axis([-r0,r0,-r0,r0])
hold off
end