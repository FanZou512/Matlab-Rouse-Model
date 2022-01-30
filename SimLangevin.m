function [R1, pos] = SimLangevin(R1,b,D,t,dt,r,Ns,constraints,N_sample)
% Langevin equation simulation, in Smoluchowskis limit
% R0 is initial position of all the beads
% A is interaction matrix
% r is friction, b is Standard Deviation of distance of neighboor beads
% ka is dn/b^2, dn is dimension
% t is time to simulate, in unit of seconds, dt is simulation time step in seconds
dimen=size(R1,2);
N=sum(Ns);
Ncons=size(R1,1)-N;
Nempty2=cumsum(Ns); % end of each chromosome
Nempty1=Nempty2-Ns+1; %start of each chromosome
Nempty1(1,:)=[];
Nempty2(end,:)=[];
%ka=dn/b^2;
%D is defined to be 1/r;
Dkt=dt*D*dimen/(b^2); % constant of D #1
Dsr2=sqrt(2*D)*sqrt(dt); % constant of D #2, 0.1 is adjust for white noise power
Det=dt*D; % constant of D #3, e is replaced by D*dt
r_min=2^(1/6)*(b/7/3); % b/7 is bead size, cutoff size is 2^1/6 larger
% initial Gaussian noise
maxkt=floor(t/dt)+1;
pos=zeros(maxkt,length(N_sample),3);
disptime=round(maxkt/40);
% core function, repeat Langevin equation, in time series
for p=0:maxkt
    kt=p*dt;
    if mod(p,disptime)==0
        disp(['Current time ',num2str(round(kt/t*100)),'%'])
    end
    w0=[wgn(N,3,0);zeros(Ncons,3)];
    f1=R1(1:N-2,:)+R1(3:N,:)-2*R1(2:N-1,:); % displacement by force of Rouse potential
    f1=[R1(2,:)-R1(1,:);f1;R1(N-1,:)-R1(N,:);zeros(Ncons,3)];
    f1(Nempty2,:)=f1(Nempty2,:)-(R1(Nempty2+1,:)-R1(Nempty2,:));
    f1(Nempty1,:)=f1(Nempty1,:)-(R1(Nempty1-1,:)-R1(Nempty1,:));
    %f_LJ=pairwiseLJ(R1,N,Det,r_min); % dU*dt/dR
    R0=R1;
    R1=R1+Dsr2*w0+Dkt*f1;%-f_LJ
    R1=PairwiseSpring(R1,Dkt,constraints);
    [R1,r]=CheckSphere(r,R0,R1); % or use CheckSphere2(r,R1);
    pos(p+1,:,1)=R1(N_sample,1)';
    pos(p+1,:,2)=R1(N_sample,2)';
    pos(p+1,:,3)=R1(N_sample,3)';
    r=max(r,1);
end
PlotRouseChain(R1,Ns,1);
end

function R=RouseSpring(R,Dkt)
% Spring potential, U(r)=kapa*(R-R0)
% index is [x-th, y-th, k], x is the index of affected bead, y is
% influencing bead, k is relative strength to Rouse
m=size(index,1);
for p=1:m
    k=Dkt*index(p,3);
    x=index(p,1);
    y=index(p,2);
    dR=(R(y,:)-R(x,:))*k;
    R(x,:)=R(x,:)+dR;
end
end

function R=PairwiseSpring(R,Dkt,index)
% Spring potential, U(r)=kapa*(R-R0)
% index is [x-th, y-th, k], x is the index of affected bead, y is
% influencing bead, k is relative strength to Rouse
m=size(index,1);
for p=1:m
    k=Dkt*index(p,3);
    x=index(p,1);
    y=index(p,2);
    dR=(R(y,:)-R(x,:))*k;
    R(x,:)=R(x,:)+dR;
end
end

function dR=pairwiseLJ(R,N,e,c)
% Lennard Jones potential, U(r)=4e[s^12/r^12-c^6/r^6], cutoff at min potential r_min=2^(1/6)*c
% Remove attraction aprt, use U(r)-U(r_min)=4e[s^12/r^12-c^6/r^6+3/4]
% derivitive: U'(r)=24e*c^6[1-2c^6/r^6]*R/r^8, e=kB*T, k1=12e*k2, k2=2c2^3,
% r2=r^2, c=r_min, r_min2=c^2
%   U'(r)=k1[1-k2/r2^3]*(R-Ri)/r2^4
% N is the max number of beads to consider repulsion [m,~]=size(R);
dR=zeros(size(R));
c2=c^2;
k2=2*c2^3;k1=12*e*k2;
for k=1:N-2
    for q=k+2:N
        x=sum((R(k,:)-R(q,:)).^2);
        if x<c2
            r2=k1*(1-k2./(x.^3))./(x.^4);
            f=r2.*(R(k,:)-R(q,:));
            if abs(f)>0.01
                f=sign(f)*0.01;
            end
            dR(k,:)=dR(k,:)+f;
            dR(q,:)=dR(q,:)-f;
        end
    end
end
end

function [R1,a]=CheckSphere2(r,R1)
rR=(sum(R1.^2,2));
b=rR>(r^2);
a=sqrt(max(rR));
R1(b,:)=R1(b,:)*(r-0.01)/a;
a=min(r,a);
end

function [R1,r]=CheckSphere(r_confine,R0,R1)
% find the first outlier after st
% BC - is a single sphere boundary
% R0 is the original pos which should be in radius, R1 is current pos
% r_confine =BC.a is the radius
r_confine_sqr=r_confine^2;
Rsq=sum((R1).^2,2);
n=(Rsq>r_confine_sqr);
c=size(R1,2);
r=r_confine;
if sum(n)>0
    p1=R0(n,:);
    p2=R1(n,:);
    p1p2=p2-p1;
    au=sum(p1p2.^2,2);
    bu=sum(p1.*p1p2,2);    % b=2*op1.*p1p2, here bu=b/2
    cu=abs(r_confine_sqr-sum(p1.^2,2));   % c=op1^2-R^2, here cu=-c
    if sum(sign(cu)==-1)>0
        disp('error on cu')
    end
    u=(sqrt(bu.^2+au.*cu)-bu)./au; % (-b+sqrt(b^2-4*a*c))/2a
    if ~isreal(sum(u(:)))
        disp('error on u')
    end
    ps=p1+u.*p1p2;       % intersection point
    xs=2*repmat(sum((p2-ps).*ps,2),1,c).*ps./repmat(sum(ps.^2,2),1,c); % shift
    R1(n,:)=R1(n,:)-xs;
    r=sqrt(max(sum((R1).^2,2)));
end
end

function [n,xs]=CheckEllipsoid(BC,x,st)
% desscription same as sphere BC
% transform into unit circle center at origin
[m,~]=size(x);
x=(x-ones(m,1)*BC.c)./(ones(m,1)*BC.a);
r=sum((x).^2,2);
% find the reflection as circle
n=find(r(st:m)>1, 1, 'first');
if isempty(n)
    xs=[];
else
    n=n(1)+st-1;
    % find the intersection
    p1=x(n-1,:);
    p2=x(n,:);
    p1p2=p2-p1;
    au=sum(p1p2.^2);
    bu=sum(p1.*p1p2);   % b=2*op1.*p1p2, here bu=b/2
    cu=1-sum(p1.^2);  % c=op1^2-R^2, here cu=-c
    u=(sqrt(bu.^2+au*cu)-bu)/au; % (-b+sqrt(b^2-4*a*c))/2a
    ps=p1+u.*p1p2;      % intersection point
    xs=-2*nansum((p2-ps).*ps).*ps/nansum(ps.^2); % shifts in transform space
    xs=xs.*BC.a;        % reverse transform
end
end
