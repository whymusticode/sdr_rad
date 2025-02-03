% clear
% function testPlanningFM
load fmTow

pixels = 4e2;
plotRange = 100; % km to plot

h0 = 0;

% lat0 = 42.460007; % lab
% lon0 = -71.2674996;
lat0 = 42.374254; % home
lon0 = -71.2476313;

Gtu = 0; % tx gain  at target [db]
Gtd = 0; % tx gain  at reciever [db]
BW = .256e6; % TV bandwidth in Hz
tau = 1; % seconds integration time
L = 5; % system losses db
RCS = 20; % dbsm
cancellation = 35; % [db] direct path cancellation
T = 1000; % noise temp K
K = 1.38e-23; % boltzmann J/K

%% end user parameters
noiseDB = 10*log10(K*T*BW); % power of noise floor

targetDB = @(Pt,r1,r2,l) 10*log10(Pt) + Gtu - L + 20*log10(l) - ...
  30*log10(4*pi) - 20*log10(r2) + RCS +10*log10(tau*BW) - ...
  20*log10(r1);

directDB = @(Pt,r,l) 10*log10(Pt) + Gtd - L + 20*log10(l) - ...
  20*log10(4*pi) - 20*log10(r);

%% narrow down to best tower choices
lat = fmTow.Lat + fmTow.itu/60 + fmTow.de/3600;
lon = -(fmTow.Long + fmTow.itu1/60 + fmTow.de1/3600); % west is -
mastHeight = fmTow.RCAM;

fmFreq = fmTow.vice;
lam = 3e2./fmFreq; % wavelength (m)
pow = fmTow.ERP*1e3; % watts

[xEast,yNorth,zUp] = geodetic2enu(lat,lon,mastHeight,...
  lat0,lon0,h0,referenceEllipsoid('wgs84'));

x = linspace(-1e3*plotRange,1e3*plotRange,pixels);y = x;
[X,Y] = meshgrid(x,y); XY = [X(:),Y(:)];

R1 = permute(sqrt(sum(XY.^2,2)),[3,2,1]); % target to reciever

XYe = [xEast,yNorth];
R3 = sqrt(sum(XYe.^2,2)); % emitter to reciever

direc = directDB(pow,R3,lam);

idx = uniquify(lam,direc); % use only loudest towers
dataMat = [XYe,lam,pow,direc,zUp];
dataM = dataMat(idx,:);

%% plot relative location of towers
% figure(81)
% plot(xEast/1e3,yNorth/1e3,'.',dataM(:,1)/1e3,dataM(:,2)/1e3,'.')
% font

R2 = sqrt(sum((permute(XY,[3,2,1]) - dataM(:,1:2)).^2,2));

signal = targetDB(dataM(:,4),R1,R2,dataM(:,3));

SIR = signal - max(dataM(:,5) - cancellation , noiseDB);

sqKM = sum(SIR > 15,3) * 4*plotRange^2/pixels^2; % square kilometers
% [3e2./dataM(:,3),sqKM]
% [3e2./dataM(:,3),dataM(:,5)]


% figure(11)
% hold on
% plot(3e2./dataM(:,3),dataM(:,5)+113,'+')

% nEm = 22;

[val,sidx] = sort(sqKM,'descend');
dataM2 = dataM(sidx,:);
disp([3e2./dataM2(:,3),sqrt(val)])



return


[~,idx,~] = intersect(3e2./dataM2(:,3),[97.3,91.7,93.3]);


figure(435)
hold off
plot(0,0,'o')
hold on
for i = 1:size(dataM2,1)
  if any(abs(FC - 3e8./dataM2(i,3))<1e-5)
    plot(dataM2(i,1)/1e3,dataM2(i,2)/1e3,'gx')
    
    text(dataM2(i,1)/1e3,dataM2(i,2)/1e3,['  ',num2str(3e2./dataM2(i,3))])
  end
end
font
axis equal
% disp('range to emitter, frequency, kilimoter coverage')
[3e2./dataM2(:,3),sqrt(sum(dataM2(:,1:2).^2,2))/1e3]
% [3e2./dataM2(:,3),dataM2(:,5) - cancellation - noiseDB]




%% end used code 

% return
% 2d matrix inverse for speed
% invMe = @(tmp)[tmp(2,2,:),-tmp(1,2,:);-tmp(2,1,:),tmp(1,1,:)]./...
%   (tmp(1,1,:).*tmp(2,2,:)-tmp(1,2,:).*tmp(2,1,:));
%
% n = length(lam);
% XYnorm = XY./R1;
% pNoRot = diag([sigD,sigC].^2);
% count = zeros(pixels^2,1);
% for i = 1:n
%
%   xms = XY - XYe(i,:);
%   R2 = sqrt(sum(xms.^2,2)); % target to reciever
%   snr = targetDB(pow(i),R1,R2,lam(i));
%   detIdx = snr > reqSNR & R1 + R2 - R3(i) > minDetRange;
%
% %   figure(19)
%   sir = SIR(R1,R2,R3(i));
%   disp(median(sir))
% %   imagesc(x/1e3,y/1e3,vec2mat(sir,pixels)');
% %   colorbar
% %   title(fmFreq(i))
% %   font
% %   disp(mean(sir))
% %   pause(.4)
%
%
% end
% return
% figure(41)
% rotRix = rotMat(atan2(.9,.3));
% pRot = rotRix'*pNoRot*rotRix;
% errorEllipse(pRot)
% axis equal


% traceP = sqrt(P(1,1,:) + P(2,2,:));
%
% measVariance = vec2mat(traceP(:),pixels)';
% countMat = vec2mat(count,pixels)';

% totalArea(j) = sum(sum(countMat>3));
% end

% doable = measVariance < 100 & countMat > 3;


%% plotting routine
%
% figure(41)
% hold off
% imagesc(x/1e3,y/1e3,countMat)
% caxis([3,n])
% hold on
% plot(xEast/1e3,yNorth/1e3,'r+')
% font
% axis equal
% colorbar
% xlabel('East (km)')
% ylabel('North (km)')
% title('Number of detections')
% set(gca,'Ydir','normal')
%
%
% countDx = countMat > 3;
% measVariance(~countDx) = 1e7;
% figure(42)
% hold off
% imagesc(x/1e3,y/1e3,measVariance)
% caxis([0,400])
% hold on
% plot(xEast/1e3,yNorth/1e3,'r+')
% font
% axis equal
% colorbar
% xlabel('East (km)')
% ylabel('North (km)')
% title('2d Measurement Accuracy (m)')
% set(gca,'Ydir','normal')
% cm = colormap;
% cm(end,:) = 1;
% colormap(cm);



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


% lat = lat(narrow);
% lon = lon(narrow);
% IDX = IDX(narrow);

% narrow = fmTow.ERP > 10;
% IDX = 1:length(fmTow.ERP);


%   for j = 1:size(measCov,3)
%     errorEllipse(measCov(:,:,j))
%     drawnow
%   end


% geometricDiversityThresh = 10; % km emitters must be this far apart
% minDetRange = 1e3; % minimum detectable range (meters)

% reqSNR = 13; % db for detection


%
% sigD = 50; % range accuracy
% sigC = 1e4; % stand in for no crossrange accuracy