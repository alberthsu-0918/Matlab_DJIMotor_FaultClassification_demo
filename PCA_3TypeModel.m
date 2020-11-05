%% new MCM model, to determine 3 cluster of air-pressor status: good, bad,
%% and before error happen

% clc;
% Dajan Test 2.1. Finding 10 peaks of wave's spectrum. 
% Apply PCA to data set of the peak and get new "projection".
% The result can be used as basics of classification.
% Refer to
% http://research.med.helsinki.fi/corefacilities/proteinchem/pca_introduction_basics.pdf
% and
% https://learnche.org/pid/latent-variable-modelling/principal-component-analysis/pca-example-analysis-of-spectral-data
% Albert Hsu  2017/01/26
% ¼Ú­C ©Ò¦VµL¼Ä

clear all;

L = 4096; %load 4096 samples of wavefile
N = 30; %number of oberservation

x1 = load('g1.csv');
x2 = load('g2.csv');
x3 = load('g3.csv');
x4 = load('g4.csv');
x5 = load('g5.csv');
x6 = load('g6.csv');
x7 = load('g7.csv');
x8 = load('g8.csv');
x9 = load('g9.csv');
x10 =load('g10.csv');


x11 = load('bad1.csv');
x12 = load('bad2.csv');
x13 = load('bad3.csv');
x14 = load('bad4.csv');
x15 = load('bad5.csv');
x16 = load('bad6.csv');
x17 = load('bad7.csv');
x18 = load('bad8.csv');
x19 = load('bad9.csv');
x20 = load('bad10.csv');

x21 = load('bfErr1.csv');
x22 = load('bfErr2.csv');
x23 = load('bfErr3.csv');
x24 = load('bfErr4.csv');
x25 = load('bfErr5.csv');
x26 = load('bfErr6.csv');
x27 = load('bfErr7.csv');
x28 = load('bfErr8.csv');
x29 = load('bfErr9.csv');
x30 = load('bfErr10.csv');

x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30];
%FrameN = 12; %length(x1)/L;
FrameN = 1;

h = hann(L);
%h = 1;
x1_ps = zeros(L/2,1); 
x2_ps = zeros(L/2,1);
x3_ps = zeros(L/2,1); 
x4_ps = zeros(L/2,1);
x5_ps = zeros(L/2,1); 
x6_ps = zeros(L/2,1);
x7_ps = zeros(L/2,1); 
x8_ps = zeros(L/2,1);
x9_ps = zeros(L/2,1); 
x10_ps = zeros(L/2,1);
x11_ps = zeros(L/2,1); 
x12_ps = zeros(L/2,1);
x13_ps = zeros(L/2,1); 
x14_ps = zeros(L/2,1);
x15_ps = zeros(L/2,1); 
x16_ps = zeros(L/2,1);
x17_ps = zeros(L/2,1); 
x18_ps = zeros(L/2,1);
x19_ps = zeros(L/2,1); 
x20_ps = zeros(L/2,1);
x21_ps = zeros(L/2,1);
x22_ps = zeros(L/2,1);
x23_ps = zeros(L/2,1);
x24_ps = zeros(L/2,1);
x25_ps = zeros(L/2,1);
x26_ps = zeros(L/2,1);
x27_ps = zeros(L/2,1);
x28_ps = zeros(L/2,1);
x29_ps = zeros(L/2,1);
x30_ps = zeros(L/2,1);

for i = 1:FrameN
   x1h = x1( (i-1)*L+1 : L*i );
   x1h = x1h.*h;
   x1_ps_temp = abs(fft(x1h));
   x1_ps_temp = x1_ps_temp(1:L/2);
   x1_ps = x1_ps + x1_ps_temp;
   %
   x2h = x2( (i-1)*L+1 : L*i );
   x2h = x2h.*h;
   x2_ps_temp = abs(fft(x2h));
   x2_ps_temp = x2_ps_temp(1:L/2);
   x2_ps = x2_ps + x2_ps_temp;
   %
   x3h = x3( (i-1)*L+1 : L*i );
   x3h = x3h.*h;
   x3_ps_temp = abs(fft(x3h));
   x3_ps_temp = x3_ps_temp(1:L/2);
   x3_ps = x3_ps + x3_ps_temp;
   %
   x4h = x4( (i-1)*L+1 : L*i );
   x4h = x4h.*h;
   x4_ps_temp = abs(fft(x4h));
   x4_ps_temp = x4_ps_temp(1:L/2);
   x4_ps = x4_ps + x4_ps_temp;
   %
   x5h = x5( (i-1)*L+1 : L*i );
   x5h = x5h.*h;
   x5_ps_temp = abs(fft(x5h));
   x5_ps_temp = x5_ps_temp(1:L/2);
   x5_ps = x5_ps + x5_ps_temp;
   %
   x6h = x6( (i-1)*L+1 : L*i );
   x6h = x6h.*h;
   x6_ps_temp = abs(fft(x6h));
   x6_ps_temp = x6_ps_temp(1:L/2);
   x6_ps = x6_ps + x6_ps_temp;
   %
   x7h = x7( (i-1)*L+1 : L*i );
   x7h = x7h.*h;
   x7_ps_temp = abs(fft(x7h));
   x7_ps_temp = x7_ps_temp(1:L/2);
   x7_ps = x7_ps + x7_ps_temp;
   %
   x8h = x8( (i-1)*L+1 : L*i );
   x8h = x8h.*h;
   x8_ps_temp = abs(fft(x8h));
   x8_ps_temp = x8_ps_temp(1:L/2);
   x8_ps = x8_ps + x8_ps_temp;
   %
   x9h = x9( (i-1)*L+1 : L*i );
   x9h = x9h.*h;
   x9_ps_temp = abs(fft(x9h));
   x9_ps_temp = x9_ps_temp(1:L/2);
   x9_ps = x9_ps + x9_ps_temp;
   %
   x10h = x10( (i-1)*L+1 : L*i );
   x10h = x10h.*h;
   x10_ps_temp = abs(fft(x10h));
   x10_ps_temp = x10_ps_temp(1:L/2);
   x10_ps = x10_ps + x10_ps_temp;
   %
   x11h = x11( (i-1)*L+1 : L*i );
   x11h = x11h.*h;
   x11_ps_temp = abs(fft(x11h));
   x11_ps_temp = x11_ps_temp(1:L/2);
   x11_ps = x11_ps + x11_ps_temp;
   %
   x12h = x12( (i-1)*L+1 : L*i );
   x12h = x12h.*h;
   x12_ps_temp = abs(fft(x12h));
   x12_ps_temp = x12_ps_temp(1:L/2);
   x12_ps = x12_ps + x12_ps_temp;
   %
   x13h = x13( (i-1)*L+1 : L*i );
   x13h = x13h.*h;
   x13_ps_temp = abs(fft(x13h));
   x13_ps_temp = x13_ps_temp(1:L/2);
   x13_ps = x13_ps + x13_ps_temp;
   %
   x14h = x14( (i-1)*L+1 : L*i );
   x14h = x14h.*h;
   x14_ps_temp = abs(fft(x14h));
   x14_ps_temp = x14_ps_temp(1:L/2);
   x14_ps = x14_ps + x14_ps_temp;
   %
   x15h = x15( (i-1)*L+1 : L*i );
   x15h = x15h.*h;
   x15_ps_temp = abs(fft(x15h));
   x15_ps_temp = x15_ps_temp(1:L/2);
   x15_ps = x15_ps + x15_ps_temp;
   %
   x16h = x16( (i-1)*L+1 : L*i );
   x16h = x16h.*h;
   x16_ps_temp = abs(fft(x16h));
   x16_ps_temp = x16_ps_temp(1:L/2);
   x16_ps = x16_ps + x16_ps_temp;
   %
   x17h = x17( (i-1)*L+1 : L*i );
   x17h = x17h.*h;
   x17_ps_temp = abs(fft(x17h));
   x17_ps_temp = x17_ps_temp(1:L/2);
   x17_ps = x17_ps + x17_ps_temp;
   %
   x18h = x18( (i-1)*L+1 : L*i );
   x18h = x18h.*h;
   x18_ps_temp = abs(fft(x18h));
   x18_ps_temp = x18_ps_temp(1:L/2);
   x18_ps = x18_ps + x18_ps_temp;
   %
   x19h = x19( (i-1)*L+1 : L*i );
   x19h = x19h.*h;
   x19_ps_temp = abs(fft(x19h));
   x19_ps_temp = x19_ps_temp(1:L/2);
   x19_ps = x19_ps + x19_ps_temp;
   %
   x20h = x20( (i-1)*L+1 : L*i );
   x20h = x20h.*h;
   x20_ps_temp = abs(fft(x20h));
   x20_ps_temp = x20_ps_temp(1:L/2);
   x20_ps = x20_ps + x20_ps_temp;
   %
   x21h = x21( (i-1)*L+1 : L*i );
   x21h = x21h.*h;
   x21_ps_temp = abs(fft(x21h));
   x21_ps_temp = x21_ps_temp(1:L/2);
   x21_ps = x21_ps + x21_ps_temp;
   %
   x22h = x22( (i-1)*L+1 : L*i );
   x22h = x22h.*h;
   x22_ps_temp = abs(fft(x22h));
   x22_ps_temp = x22_ps_temp(1:L/2);
   x22_ps = x22_ps + x22_ps_temp;
   %
   x23h = x23( (i-1)*L+1 : L*i );
   x23h = x23h.*h;
   x23_ps_temp = abs(fft(x23h));
   x23_ps_temp = x23_ps_temp(1:L/2);
   x23_ps = x23_ps + x23_ps_temp;
   %
   x24h = x24( (i-1)*L+1 : L*i );
   x24h = x24h.*h;
   x24_ps_temp = abs(fft(x24h));
   x24_ps_temp = x24_ps_temp(1:L/2);
   x24_ps = x24_ps + x24_ps_temp;
   %
   x25h = x25( (i-1)*L+1 : L*i );
   x25h = x25h.*h;
   x25_ps_temp = abs(fft(x25h));
   x25_ps_temp = x25_ps_temp(1:L/2);
   x25_ps = x25_ps + x25_ps_temp;
   %
   x26h = x26( (i-1)*L+1 : L*i );
   x26h = x26h.*h;
   x26_ps_temp = abs(fft(x26h));
   x26_ps_temp = x26_ps_temp(1:L/2);
   x26_ps = x26_ps + x26_ps_temp;
   %
   x27h = x27( (i-1)*L+1 : L*i );
   x27h = x27h.*h;
   x27_ps_temp = abs(fft(x27h));
   x27_ps_temp = x27_ps_temp(1:L/2);
   x27_ps = x27_ps + x27_ps_temp;
   %
   x28h = x28( (i-1)*L+1 : L*i );
   x28h = x28h.*h;
   x28_ps_temp = abs(fft(x28h));
   x28_ps_temp = x28_ps_temp(1:L/2);
   x28_ps = x28_ps + x28_ps_temp;
   %
   x29h = x29( (i-1)*L+1 : L*i );
   x29h = x29h.*h;
   x29_ps_temp = abs(fft(x29h));
   x29_ps_temp = x29_ps_temp(1:L/2);
   x29_ps = x29_ps + x29_ps_temp;
   %
   x30h = x30( (i-1)*L+1 : L*i );
   x30h = x30h.*h;
   x30_ps_temp = abs(fft(x30h));
   x30_ps_temp = x30_ps_temp(1:L/2);
   x30_ps = x30_ps + x30_ps_temp;
   %
end

%normalize and get rid off offset of power spectrum
%x_ps = [x1_ps-mean(x1_ps) x2_ps-mean(x2_ps) x3_ps-mean(x3_ps) x4_ps-mean(x4_ps) x5_ps-mean(x5_ps) x6_ps-mean(x6_ps) x7_ps-mean(x7_ps) x8_ps-mean(x8_ps) x9_ps-mean(x9_ps) x10_ps-mean(x10_ps) x11_ps-mean(x11_ps) x12_ps-mean(x12_ps) x13_ps-mean(x13_ps) x14_ps-mean(x14_ps) x15_ps-mean(x5_ps) x16_ps-mean(x16_ps) x17_ps-mean(x17_ps) x18_ps-mean(x18_ps) x19_ps-mean(x19_ps) x20_ps-mean(x20_ps)];
x_ps = [x1_ps x2_ps x3_ps x4_ps x5_ps x6_ps x7_ps x8_ps x9_ps x10_ps x11_ps x12_ps x13_ps x14_ps x15_ps x16_ps x17_ps x18_ps x19_ps x20_ps x21_ps x22_ps x23_ps x24_ps x25_ps x26_ps x27_ps x28_ps x29_ps x30_ps];
% x_ps = [x1_ps/sum(x1_ps) x2_ps/sum(x2_ps) x3_ps/sum(x3_ps) x4_ps/sum(x4_ps) x5_ps/sum(x5_ps) x6_ps/sum(x6_ps) x7_ps/sum(x7_ps) x8_ps/sum(x8_ps) x9_ps/sum(x9_ps) x10_ps/sum(x10_ps) x11_ps/sum(x11_ps) x12_ps/sum(x12_ps) x13_ps/sum(x13_ps) x14_ps/sum(x14_ps) x15_ps/sum(x15_ps) x16_ps/sum(x16_ps) x17_ps/sum(x17_ps) x18_ps/sum(x18_ps) x19_ps/sum(x19_ps) x20_ps/sum(x20_ps)];
%x_ps = [x1_ps/sum(x1_ps) x2_ps/sum(x2_ps) x3_ps/sum(x3_ps) x4_ps/sum(x4_ps) x5_ps/sum(x5_ps) x6_ps/sum(x6_ps) x7_ps/sum(x7_ps) x8_ps/sum(x8_ps) x9_ps/sum(x9_ps) x10_ps/sum(x10_ps)];

x_ps = 20*log10(x_ps);
scaleMean = zeros(L/2,1);
for i = 1:L/2
    tempXs = x_ps(i,:);
    scaleMean(i) = mean(tempXs);
    x_ps(i,:) = tempXs - mean(tempXs);
end

% for i = 1:N
%     tempXs = x_ps(:,i);
%     x_ps(:,i) = tempXs - mean(tempXs);
% end


%% Find 10 peaks as characteristics 
%% useless
% psSum = zeros(L/2,1);
% for i = 1:N
%     psSum = psSum + x_ps(:,N);
% end
% 
% loc = peakfinder(psSum);
% x_ps_peak = [x_ps(loc,1) x_ps(loc,2) x_ps(loc,3) x_ps(loc,4) x_ps(loc,5) x_ps(loc,6) x_ps(loc,7) x_ps(loc,8) x_ps(loc,9) x_ps(loc,10) ];


%PCA Analysis
[pcs, trans, evs] = princomp(x_ps');
%[pcs, trans, evs] = princomp(x_ps_peak');
%[trans,pcs,evs] = pca1(x_ps_peak);
% trans = real(trans);
% Below are SVD transpose
pcsN = 10;
[U , S, V] = svd(x_ps');
%trans = (pcs'*x_ps)';
trans = (V'*x_ps)';


plot(trans(1,1:pcsN))
hold on
plot(trans(2,1:pcsN))
plot(trans(3,1:pcsN))
plot(trans(4,1:pcsN))
plot(trans(5,1:pcsN))
plot(trans(6,1:pcsN))
plot(trans(7,1:pcsN))
plot(trans(8,1:pcsN))
plot(trans(9,1:pcsN))
plot(trans(10,1:pcsN))
plot(trans(11,1:pcsN), 'r')
plot(trans(12,1:pcsN), 'r')
plot(trans(13,1:pcsN), 'r')
plot(trans(14,1:pcsN), 'r')
plot(trans(15,1:pcsN), 'r')
plot(trans(16,1:pcsN), 'r')
plot(trans(17,1:pcsN), 'r')
plot(trans(18,1:pcsN), 'r')
plot(trans(19,1:pcsN), 'r')
plot(trans(20,1:pcsN), 'r')
plot(trans(21,1:pcsN), 'g')
plot(trans(22,1:pcsN), 'g')
plot(trans(23,1:pcsN), 'g')
plot(trans(24,1:pcsN), 'g')
plot(trans(25,1:pcsN), 'g')
plot(trans(26,1:pcsN), 'g')
plot(trans(27,1:pcsN), 'g')
plot(trans(28,1:pcsN), 'g')
plot(trans(29,1:pcsN), 'g')
plot(trans(30,1:pcsN), 'g')

score = [trans(1,1:2); trans(2,1:2) ;trans(3,1:2); trans(4,1:2); trans(5,1:2); trans(6,1:2); trans(7,1:2); trans(8,1:2); trans(9,1:2); trans(10,1:2)];

% t = 1:51200/L:51200/2;

% %% SVD Test
% for i = 1:5
%     meanV(:,i) = mean(x(:,i))*ones(L,1);  
% end
% xz = x - meanV;
% %[U , S, V] = svds(xz,SensorN);
% [U , S, V] = svds(xz,1);
% xr = U*S*V';
% xr_fs = abs(fft(xr(:,1)));
% xr_fs = xr_fs(1:length(xr_fs)/2);



