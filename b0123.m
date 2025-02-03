clear
fc = 97.9e6;
fs = 1e6;
g = 30;
N = 1e7;
dataType = 'int16';
save('params')


delete 'data\*.mat'
!matlab -r b0 &
!matlab -r b1 &
!matlab -r b2 &
!matlab -r b3 &