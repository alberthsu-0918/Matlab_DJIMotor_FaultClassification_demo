% clc;
% MCM Test for compressor 
% Apply PCA to data set of the peak and get new "projection".
% The result can be used as basics of classification.
% Refer to
% http://research.med.helsinki.fi/corefacilities/proteinchem/pca_introduction_basics.pdf
% and
% https://learnche.org/pid/latent-variable-modelling/principal-component-analysis/pca-example-analysis-of-spectral-data
% Albert Hsu  2017/11/23
% ¼Ú­C ©Ò¦VµL¼Ä

clear all;

L = 4096; %load 4096 samples of wavefile
N = 92; %number of oberservation
FrameN = 5;  %for sound data, FrameN=5

OKDir = 'C:\Matlab\work\PCA\MCM\sound2017_05_17\';
%OKDir = 'C:\Matlab\work\PCA\MCM\sound2017_05_17\';
 %OKDir = 'C:\Matlab\work\PCA\MCM\2017_05_08to09\';
% OKDir = 'C:\Matlab\work\PCA\MCM\2017_11_22_afRepair\'; 
% OKDir = 'C:\Matlab\work\PCA\MCM\2017_0615_afRepair\';

temp0=dir([OKDir,'*.csv']);
num_temp0=length(temp0);

failDir = 'C:\Matlab\work\PCA\MCM\sound2017_06_03\';
%failDir = 'C:\Matlab\work\PCA\MCM\2017_0603_bfFail\'; %2017_0603_bfFail\one day bf 1st breaking day     %2017_11_20\:one day bf 2nd breaking day   %2017_06_04_TheBreakdingDAy
temp1=dir([failDir,'*.csv']);
num_temp1=length(temp1);

N = num_temp0 + num_temp1; 
x = zeros(L*FrameN , N);


for i = 1:num_temp0
    filename=[OKDir,temp0(i).name];
    x(:, i) = load(filename);
end

for i = 1:num_temp1
    filename=[failDir,temp1(i).name];
    x(:, num_temp0+i) = load(filename);
end




%FrameN = 12; %length(x1)/L;


h = hann(L);

x_ps = zeros(L/2,N);
x_h =  zeros(L,N);


% average FFT spectrum
for j = 1:N
    for i = 1:FrameN
        x_h(:,j) = x((i-1)*L+1 : L*i, j);
        x_h(:,j) = x_h(:,j).*h;
        x_ps_temp = abs(fft(x_h(:,j)));
        x_ps_temp = x_ps_temp(1:L/2);
        x_ps(:,j) = x_ps(:,j) + x_ps_temp; 
    end
end



%normalize and get rid off offset of power spectrum
x_ps = 20*log10(x_ps);
x_pstest = x_ps;
scaleMean = zeros(L/2,1);
for i = 1:L/2
    tempXs = x_ps(i,:);
    scaleMean(i) = mean(tempXs);
    x_ps(i,:) = tempXs - mean(tempXs);
end

%normalized std for each sample
% for i = 1:N
%     x_ps(:,i) = x_ps(:,i)./ sqrt(sum(x_ps(:,i).*x_ps(:,i))); 
%     std = sqrt(sum(x_ps(:,i).*x_ps(:,i))); 
% end

% normalized std for each spectrum data
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./ sqrt(sum(x_ps(i,:).*x_ps(i,:))); 
end

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
%[pcs, trans, evs] = princomp(x_ps');
%[pcs, trans, evs] = princomp(x_ps_peak');
%[trans,pcs,evs] = pca1(x_ps_peak);
% trans = real(trans);
% Below are SVD transpose

pcsN = 10;
[U , S, V] = svd(x_ps', 0);

%trans = (pcs'*x_ps)';
trans = (V'*x_ps)';


OKnum = num_temp0;
failnum = num_temp1;

plot(trans(1,1:pcsN))
hold on
for i = 2:OKnum
    plot(trans(i,1:pcsN))
end

for i = OKnum+1 : failnum + OKnum 
    plot(trans(i,1:pcsN), 'r');
end


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



