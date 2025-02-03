
% load from two seperate runs, data and DATA (2 different channels) 

dat = [data,DATA];

tmp = [];
for i = 1:size(dat,3)
  tmp = [tmp;dat(:,:,i)];
end
DATA = tmp;

data = DATA(1:7e5,1);
data2 = DATA(1:7e5,2);
[val,lag] = xcorr(data,data2);
figure(45)
plot(lag,abs(val))

[~,loc] = max(abs(val)); LAG = lag(loc);
if LAG > 0;D1 = DATA(LAG:end,1); D2 = DATA(1:end-LAG+1,2);
else;LAG = -LAG; D2 = DATA(LAG:end,2); D1 = DATA(1:end-LAG+1,1);
end


data = D1(1:7e5);
data2 = D2(1:7e5);
data = D1(end-7e5:end);
data2 = D1(1:7e5);
[val,lag] = xcorr(data,data2);
figure(45)
plot(lag,abs(val))

DATA = [D1,D2];

clear D1 D2 tmp lag val dat data data2
save('5_26_374Secs')


% data = DATA(end-7e5:end,1);
% data2 = DATA(end-7e5:end,2);
data = DATA(1:7e5,1);
data2 = DATA(1:7e5,2);
[val,lag] = xcorr(data,data2);
figure(45)
plot(lag,abs(val))

