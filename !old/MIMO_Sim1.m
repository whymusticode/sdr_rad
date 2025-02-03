%MIMO_Sim1.m
%Writen by Brad Herrick 3/15/18
%Only 2D version used here
%Considers multiple radars operating in a MIMO fassion collecting all
%bistatic ranges
%Collects number of intersections of ellipses for detection (Hit)
%Animation option by setting Nt>1;
%boldify_17.m will be used if the file is available
% The 1980's tried and true boldify.m by Steve Smith does not work in the
% latest matlab.  Either get an update or use the boldify_17 included here
% with the update

clear
% close all
%% Waveform
WF.FGHz=.5;
WF.wl=.3./WF.FGHz;
WF.PRF=500;
WF.Vamb=WF.wl.*WF.PRF./2;
WF.Tcpi=.5;
WF.Npulse=floor(WF.Tcpi.*WF.PRF);
WF.NCPI=2;  %Length of simulated data
WF.Nfft=16.*512;
WF.Vres=WF.Vamb./WF.Nfft;
WF.Rres=30;

% Set of Radars
Radar.N=5;
Radar.Rad=1000;
Radar.Phi=2.*pi.*[1:Radar.N]./Radar.N;
Radar.x0=1000+Radar.Rad.*cos(Radar.Phi);
Radar.y0=1000+Radar.Rad.*sin(Radar.Phi);
%Radar.z0=zeros(size(Radar.x0));

%% Time steps
dt=10;
Nt=200;
t=[0:Nt-1].*dt;
%% Simulate scenario and absolute positions
Radar.vx=ones(1,Radar.N).*0;
Radar.vy=ones(1,Radar.N).*0;
%Radar.vz=ones(1,Radar.N).*0;
Target.N=10;
Target.xc=3000;
Target.yc=3000;
Target.sigmax=500;
Target.sigmay=500;
Target.case='random';
switch Target.case
  case 'circle'
    Theta=[1:Target.N]./Target.N.*2.*pi;
    Target.x0=Target.xc+Target.sigmax.*cos(Theta);
    Target.y0=Target.yc+Target.sigmay.*sin(Theta);
  case 'random'
    Target.x0=Target.xc+Target.sigmax.*randn(1,Target.N);
    Target.y0=Target.xc+Target.sigmay.*randn(1,Target.N);
end
%Window: pixels are set to 1/2 range resolution
xs=[0:WF.Rres./2:4000];
ys=[0:WF.Rres./2:4000];
[Xs,Ys]=meshgrid(xs,ys);
Target.X=Target.x0;
Target.Y=Target.y0;

%Simulation
for it=1:Nt
  %disp(['t = ',num2str(t(it)),' (sec)'])
  Radar.X=Radar.x0+t(it)*Radar.vx;
  Radar.Y=Radar.y0+t(it)*Radar.vy;
  %Radar.z=Radar.z0+t.*Radar.vz;
  Hit=zeros(size(Xs));
  %Random walk
  Target.vx=2.*randn(1,Target.N)-1;
  Target.vy=2.*randn(1,Target.N)-1;
  for k=1:Target.N
    Target.X(k)=Target.X(k)+dt*Target.vx(k);
    Target.Y(k)=Target.Y(k)+dt*Target.vy(k);
    R0=sqrt((Target.X(k)-Radar.X).^2+(Target.Y(k)-Radar.Y).^2);
    [R1,R2]=meshgrid(R0,R0);
    RBS(:,:)=R1+R2;   %Bistatic Range;
    
    %% Brute force: Look for overlaps in bistatic range
    for i=1:Radar.N
      Ri=sqrt((Xs-Radar.X(i)).^2+(Ys-Radar.Y(i)).^2);
      for j=1:Radar.N
        Rj=sqrt((Xs-Radar.X(j)).^2+(Ys-Radar.Y(j)).^2);
        Rij=Ri+Rj;
        ihit=find(abs(Rij-RBS(i,j)) <= WF.Rres);
        Hit(ihit)=Hit(ihit)+1;
      end
    end
  end
  N2=Radar.N.^2;
  isbad=find(Hit<N2-Radar.N);
  HitDet=Hit;
  HitDet=min(HitDet,N2+Radar.N);
  HitDet(isbad)=nan;
  %Detect=Cluster(Hit,Thresh);
  
  %% Plotting
%   set(gcf,'position',[1048 431 858 680])
  pcolor(xs./1000,ys./1000,HitDet)
  colormap jet
  shading flat
  axis([0 4 0 4])
  colorbar
  axis xy
  hold on
  plot(Target.X./1000,Target.Y./1000,'or')
  plot(Radar.X./1000,Radar.y0./1000,'+k')
  xlabel('X (km)')
  ylabel('Y (km)')
  title(['Overlap of Bistatic Ranges, ',num2str(WF.Rres),...
    ' (m) Range Resolution']);
  if exist('boldify_17','file')==2
               boldify_17
  end
  hold off
  drawnow
  M(it) = getframe;
end

v = VideoWriter('multistatic.mp4','MPEG-4');
open(v);
writeVideo(v,M);
close(v);

