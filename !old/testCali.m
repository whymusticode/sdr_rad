







n = 1e4;

filter = ifft(   fft(rxa(1:n))./(fft(txa(1:n)) + 1*exp(rand(n,1)*1i*2*pi)    )    );

TX = fft(txa(1:n2)).*fft(filter,n2);
tx = ifft(TX);



n2 = 1e5;
rxs = rxa(1:n2);
RX = fft(rxs);
TX = fft(txa(1:n2));
TXs = 20*log10(smpl(abs(fftshift(TX)),1e2));
RXs = 20*log10(smpl(abs(fftshift(RX)),1e2));

x = 1:1e3;
plot(x,TXs,x,RXs)





freqIn = 2.4e6; % sample freq in
freqOut = .25e6; % band pass filter output
txn = extractChannels(tx,0,freqIn,freqOut);
rxn = extractChannels(rxs,0,freqIn,freqOut);


rr = gpuArray(single(rxa(1:3e5)));
tt = gpuArray(single(txa(1:3e5)));
rxClean = removeDirectVect(300,20,rr,tt);
10*log10(abs((rxClean'*rxClean) / (rr'*rr)))



cfs = (-10:2:10)/24;
bpFreq = .15;
txchannels = extractChannels(txa(1:2*p.fs),cfs,bpFreq);
rxchannels = extractChannels(rxa(1:2*p.fs),cfs,bpFreq);

t1 = txchannels(:,6);
r1 = rxchannels(:,6);


r1fr = smpl(log(abs(fftshift(fft(r1(1:8:end))))),100);
t1fr = smpl(log(abs(fftshift(fft(t1(1:8:end))))),100);


t1c = chanCali(t1,r1(1:1e1));

t1cfr = smpl(log(abs(fftshift(fft(t1c(1:8:end))))),100);

x = 1:length(r1fr);
plot(x,r1fr,x,t1fr,x,t1cfr)


rxClean = removeDirectVect(300,20,rr,tt);
10*log10(abs((rxClean'*rxClean) / (rr'*rr)))























