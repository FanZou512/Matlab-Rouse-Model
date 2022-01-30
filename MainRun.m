% initialize positions
Ns=[80;80;80;80;80;80]; % number of monomers in each segements
b=0.2; % average RMSD of adjacent monomer distance, unit um
D=1e-3; % diffusion constant, um^2/s
[R0,N,r0] = InitialChain(Ns,b/2); % R0 initial pos, N total number of monomers, r0 radius
t0=b^2/3/pi^2/D; % monomer relaxation, 2.5*tN
tN=b^2/12/D; % tau_N for monomer relaxation
tR=Ns(1)^2*b^2/3/pi^2/D; % tau_Rouse for chain relaxation
slope=2*b*sqrt(3*D/pi);

%% extra constraints
R_ref=[0,1,0;0,-1,0]; % tethering centromere/telomere
x=cumsum([0;Ns]);x=x(1:end-1)+1;
N_seg=size(Ns,1);
constraints=[ones(N_seg,1).*x,ones(N_seg,1)*(N+1),ones(N_seg,1)*100];
constraints=[constraints;ones(N_seg,1).*(x+Ns-1),ones(N_seg,1)*(N+2),ones(N_seg,1)*50];
R0=[R0;R_ref];
N_sample=[20;60;100;140;180;220;260;300;340;380;420;460];%N_sample=x+round(Ns/2)-1;
    %[2;40;82;120;162;200;242;280;322;360;402;440];

%% Run the simulation for relaxation
dt=0.1;
t=36000;
[R1, pos]= SimLangevin(R0,b,D,t,dt,r0+b/10,Ns,constraints,N_sample);
rmax=sqrt(max(sum(R1.^2,2)));

%% extra dynamics
t=72000;
[R1,pos ]= SimLangevin(R0,b,D,t,dt,rmax+b/10,Ns,constraints,N_sample);
rmax=sqrt(max(sum(R1.^2,2)));
PlotRouseChain(R1,Ns,1)

%% post analysis
%t=0:length(pos)-1;
%figure(),plot(t,pos)
T=unique(round((100:20:10000).^1.0));
MSD1 =MSDspecify(T,pos(:,1:2:12,1),pos(:,1:2:12,2),pos(:,1:2:12,3));
MSD2 =MSDspecify(T,pos(:,2:2:12,1),pos(:,2:2:12,2),pos(:,2:2:12,3));
MSDs=MSD1+MSD2;
RMSD =RMSDspecify(T,pos(:,:,1),pos(:,:,2),pos(:,:,3));
T=T*dt;
%figure,loglog(t,MSD)