function PlotRouseChain(R0,Ns,option)
if nargin<3
    option=1;
end
figure(),hold on
s1=1;
colorlist={'#D95319','#EDB120';'#7E2F8E','#0072BD';'#77AC30','#4DBEEE';'#A2142F','#FD00FF';...
    '#1FBF2F','#BFEFBF';'#EDB120','#D95319';'#7E2F8E','#0072BD';'#4DBEEE','#77AC30'};
for k=1:length(Ns)
    s2=s1+Ns(k)-1;
    plot(R0(s1:s2,1),R0(s1:s2,2),'-o','LineWidth',3,'Color',colorlist{k,1},'MarkerSize',...
        12,'MarkerEdgeColor',colorlist{k,1},'MarkerFaceColor',colorlist{k,2})
    s1=s1+Ns(k);
end
rmax=sqrt(max(sum(R0.^2,2)));
if option==1
    x=(-1:0.01:1);
    y=sqrt(1-x.^2);
    x=[x,x(end:-1:1)];y=[y,-y];
    x=x*rmax;y=y*rmax;
    plot(x,y,'-','LineWidth',3,'MarkerSize',12,'Color','#77AC30')
    axis([-rmax,rmax,-rmax,rmax])
end
end

