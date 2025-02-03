clear

load fmTow

pixels = 1e2;
plotRange = 100; % km to plot

sigD = 50; % range accuracy
sigC = 1e4; % stand in for no crossrange accuracy

lat0 = 42.3684594; % observer location
lon0 = -71.2459723;

BW = .256e6; % TV bandwidth in Hz 

Gtu = 0; % gain up (towards target) db
Gtd = 0; % gain down (towards reciever) db
tau = .5; % seconds integration time
L = 5; % system losses db 
RCS = 10; % dbsm

geometricDiversityThresh = 10; % km emitters must be this far apart
minDetRange = 1e3; % minimum detectable range (meters)

reqSNR = 13; % db for detection

T = 500; % noise temp K
K = 1.38e-23; % boltzmann J/K


%% end user parameters

targetSNR = @(Pt,r1,r2,l) 10*log10(Pt) + Gtu - L + 20*log10(l) - ...
  30*log10(4*pi) - 20*log10(r2) + RCS - 10*log10(K*T)+10*log10(tau) - ...
  20*log10(r1);

directSNR = @(Pt,r,l) 10*log10(Pt) + Gtd - L + 20*log10(l) - ...
  20*log10(4*pi) - 20*log10(r)-10*log10(K*T)-10*log10(BW);


SIR = @(r1,r2,r3) Gtu - Gtd  - 10*log10(4*pi) - 20*log10(r2) ...
  - 20*log10(r1) + 20*log10(r3) + RCS + 10*log10(tau) + 10*log10(BW);
% ^ with signal processing gain

%% narrow down to best tower choices
narrow = fmTow.ERP > 10;
IDX = 1:length(fmTow.ERP);
lat = fmTow.Lat + fmTow.itu/60 + fmTow.de/3600;
lon = fmTow.Long + fmTow.itu1/60 + fmTow.de1/3600;
mastHeight = fmTow.RCAM;
lat = lat(narrow);
lon = lon(narrow);
IDX = IDX(narrow);
mastHeight = mastHeight(narrow);
[erp,idx] = sort(fmTow.ERP(narrow),'descend');
IDX = IDX(idx);
lat = lat(idx);
lon = -lon(idx);
mastHeight = mastHeight(idx);
r0 = 6371; % earth radius (km)
n = length(lat);
tooClose = 0;
thresh = (geometricDiversityThresh*180/pi/r0)^2;
for i = 1:n
  tci = sum(([lat(i+1:end),lon(i+1:end)] - [lat(i),lon(i)]).^2,2) < thresh;
  tmp = [zeros(i,1);tci];
  tooClose = tooClose + tmp;
end
lat = lat(tooClose==0);
lon = lon(tooClose==0);
IDX = IDX(tooClose==0);
mastHeight = mastHeight(tooClose==0);

fmFreq = fmTow.vice(IDX);
lam = 3e2./fmFreq; % wavelength (m)
if length(unique(lam)) ~= length(lam); error('overused freq');end
pow = fmTow.ERP(IDX)*1e3; % watts 

[xEast,yNorth,zUp] = geodetic2enu(lat,lon,mastHeight,...
  lat0,lon0,0,referenceEllipsoid('wgs84'));


% plot(xEast/1e3,yNorth/1e3,'.')
% font

x = linspace(-1e3*plotRange,1e3*plotRange,pixels);
y = x;

[X,Y] = meshgrid(x,y);
XY = [X(:),Y(:)];

R1 = sqrt(sum(XY.^2,2)); % target to reciever

XYe = [xEast,yNorth];
R3 = sqrt(sum(XYe.^2,2)); % emitter to reciever 

% isolation calculations 
direc = directSNR(pow,R3,lam);
% direc - 10*log10(6e6*.2);
% dynRangReq = log2( 10.^(direc/10)/(6e6*.2) ); % -processing gain

P = repmat(1e6*eye(2),[1,1,size(XY,1)]); % initialize covariance

% 2d matrix inverse for speed
invMe = @(tmp)[tmp(2,2,:),-tmp(1,2,:);-tmp(2,1,:),tmp(1,1,:)]./...
  (tmp(1,1,:).*tmp(2,2,:)-tmp(1,2,:).*tmp(2,1,:));

n = length(lam);
XYnorm = XY./R1;
pNoRot = diag([sigD,sigC].^2);
count = zeros(pixels^2,1);
for i = 1:n
  
  xms = XY - XYe(i,:);
  R2 = sqrt(sum(xms.^2,2)); % target to reciever
  snr = targetSNR(pow(i),R1,R2,lam(i));
  detIdx = snr > reqSNR & R1 + R2 - R3(i) > minDetRange;
  
%   figure(19)
  sir = SIR(R1,R2,R3(i));
  disp(median(sir))
%   imagesc(x/1e3,y/1e3,vec2mat(sir,pixels)');
%   colorbar
%   title(fmFreq(i))
%   font
%   disp(mean(sir))
%   pause(.4)
  
  count = count + detIdx;
  
  xmsNorm = xms(detIdx,:)./R2(detIdx);
  
  rotVal = permute(xmsNorm + XYnorm(detIdx,:),[2,3,1]); % rotation vector
  rotRix = rotMat(atan2(rotVal(2,:,:),rotVal(1,:,:)));
  rotRixP = permute(rotRix,[2,1,3]); % transpose of rotation matrix
  measCov = mm3d(mm3d(rotRixP,pNoRot),rotRix);
  
  P(:,:,detIdx) = invMe(invMe(P(:,:,detIdx))+invMe(measCov));
  
end
return
% figure(41)
% rotRix = rotMat(atan2(.9,.3));
% pRot = rotRix'*pNoRot*rotRix;
% errorEllipse(pRot)
% axis equal


traceP = sqrt(P(1,1,:) + P(2,2,:));

measVariance = vec2mat(traceP(:),pixels)';
countMat = vec2mat(count,pixels)';

% totalArea(j) = sum(sum(countMat>3));
% end

% doable = measVariance < 100 & countMat > 3;


%% plotting routine

figure(41)
hold off
imagesc(x/1e3,y/1e3,countMat)
caxis([3,n])
hold on
plot(xEast/1e3,yNorth/1e3,'r+')
font
axis equal
colorbar
xlabel('East (km)')
ylabel('North (km)')
title('Number of detections')
set(gca,'Ydir','normal')


countDx = countMat > 3;
measVariance(~countDx) = 1e7;
figure(42)
hold off
imagesc(x/1e3,y/1e3,measVariance)
caxis([0,400])
hold on
plot(xEast/1e3,yNorth/1e3,'r+')
font
axis equal
colorbar
xlabel('East (km)')
ylabel('North (km)')
title('2d Measurement Accuracy (m)')
set(gca,'Ydir','normal')
cm = colormap;
cm(end,:) = 1;
colormap(cm);
% figure(43)
% hold off
% imagesc(x/1e3,y/1e3,doable)
% % hold on
% % plot(xEast/1e3,yNorth/1e3,'w+')
% font
% axis equal
% colorbar
% xlabel('East (km)')
% ylabel('North (km)')
% title('Can we expect to find a target?')
% colormap(bone)
% set(gca,'Ydir','normal')

% figure(44)
% plot(potFREQvec,totalArea/1e6*(x(2)-x(1))^2)
% font
% xlabel('Center Fequency (MHz)')
% ylabel('Coverage (km^2)')
%% end of useable code


%   for j = 1:size(measCov,3)
%     errorEllipse(measCov(:,:,j))
%     drawnow
%   end



