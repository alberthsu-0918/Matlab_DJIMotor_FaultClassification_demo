% clc;
% MCM Test for compressor.  
% After run 'wavelet_spectrum_PCAModel.m' to get pcs, the transpose matrix
% Run this m file to see if the pcs can project other raw data to a better
% 3d axis.
% Further clustering is needed
% 2018/04/02
%clear all;


L = 4096; %load 4096 samples of wavefile
%N = 20; %number of oberservation
N = 0;
FrameN = 1;
num_temp0 = 0;
num_temp1 = 0;
num_temp2 = 0;
num_temp3 = 0;
num_temp4 = 0;
features = 4;
dwtDepth = 3;
base = 'coif4';

%% 1. --------------loading Data-----------------------%%
Dir0 = 'C:\Matlab\work\PCA\MCM\2017_11_20\';
%Dir0 = 'C:\Matlab\work\PCA\MCM\sound2017_08_15_all\';
temp0=dir([Dir0,'*.csv']);
num_temp0=length(temp0);
%
% %Dir1 = 'C:\Matlab\work\PCA\MCM\2017_10_03\'; %2017_0603_bfFail\ %2017_06_04_TheBreakdingDAy\
Dir1 = 'C:\Matlab\work\PCA\MCM\2017_10_17\';
temp1=dir([Dir1,'*.csv']);
num_temp1=length(temp1);
% % % % 
Dir2 = 'C:\Matlab\work\PCA\MCM\2017_08_14to19\'; %2017_08_14to19
temp2=dir([Dir2,'*.csv']);
num_temp2=length(temp2);
% % % % 
Dir3 = 'C:\Matlab\work\PCA\MCM\2017_07_03\'; %2017_07_03
temp3=dir([Dir3,'*.csv']);
num_temp3=length(temp3);
% %
Dir4 = 'C:\Matlab\work\PCA\MCM\2017_06_20\'; %2017_0603_bfFail\ %2017_06_04_TheBreakdingDAy\
temp4=dir([Dir4,'*.csv']);
num_temp4=length(temp4);

N = num_temp0 + num_temp1 + num_temp2 + num_temp3 + num_temp4; 
x = zeros(L*FrameN,N);

for i = 1:num_temp0
    filename=[Dir0,temp0(i).name];
    x(:, i) = load(filename);
end

for i = 1:num_temp1
    filename=[Dir1,temp1(i).name];
    x(:, num_temp0+i) = load(filename);
end

for i = 1:num_temp2
    filename=[Dir2,temp2(i).name];
    x(:, num_temp0+ num_temp1+i) = load(filename);
end

for i = 1:num_temp3
    filename=[Dir3,temp3(i).name];
    x(:, num_temp0+ num_temp1+ num_temp2+i) = load(filename);
end

for i = 1:num_temp4
    filename=[Dir4,temp4(i).name];
    x(:, num_temp0+ num_temp1+ num_temp2+ num_temp3 + i) = load(filename);
end


% 2. --------------using d-wavelt transform to seperate data as 3 level-------------- 
% first level decomposition A D
[A,D] = dwt(x(:,1), base);
len1a = length(A);
len1d = length(D);
xA = zeros(len1a,N);
xD = zeros(len1d,N);
for i = 1:N
    [xA(:,i),xD(:,i)] = dwt(x(:,i), base);
end

% second level decomposition A ->AA AD , D-> DA DD
%  AA DA
[AA, DA] = dwt(xA(:,1), base);
len2aa = length(AA);
len2da = length(DA);
xAA = zeros(len2aa, N);
xDA = zeros(len2da, N);
for i = 1:N
    [xAA(:,i), xDA(:,i)] = dwt(xA(:,i), base);
end
%  AD DD
[AD, DD] = dwt(xD(:,1), base);
len2ad = length(AD);
len2dd = length(DD);
xAD = zeros(len2ad, N);
xDD = zeros(len2dd, N);
for i = 1:N
    [xAD(:,i), xDD(:,i)] = dwt(xD(:,i), base);
end

%third level decomposition AA-> AAA DAA, DA-> ADA DDA, AD-> AAD DAD, DD-> ADD DDD
%  AA-> AAA DAA
[AAA, DAA] = dwt(xAA(:,1), base);
len2aaa = length(AAA);
len2daa = length(DAA);
xAAA = zeros(len2aaa, N);
xDAA = zeros(len2daa, N);
for i = 1:N
    [xAAA(:,i), xDAA(:,i)] = dwt(xAA(:,i), base);
end

%  DA-> ADA DDA
[ADA, DDA] = dwt(xDA(:,1), base);
len2ada = length(ADA);
len2dda = length(DDA);
xADA = zeros(len2ada, N);
xDDA = zeros(len2dda, N);
for i = 1:N
    [xADA(:,i), xDDA(:,i)] = dwt(xDA(:,i), base);
end

%  AD-> AAD DAD  
[AAD, DAD] = dwt(xAD(:,1), base);
len2aad = length(AAD);
len2dad = length(DAD);
xAAD = zeros(len2aad, N);
xDAD = zeros(len2dad, N);
for i = 1:N
    [xAAD(:,i), xDAD(:,i)] = dwt(xAD(:,i), base);
end

%  DD-> ADD DDD
[ADD, DDD] = dwt(xDD(:,1), base);
len2add = length(ADD);
len2ddd = length(DDD);
xADD = zeros(len2add, N);
xDDD = zeros(len2ddd, N);
for i = 1:N
    [xADD(:,i), xDDD(:,i)] = dwt(xDD(:,i), base);
end
%%%%%%%%%%% end of d-wavelt transform %%%%%%%%%%%


%%use AAA as new raw data
L = length(xAAA);
h = hann(L);

x_ps = zeros(L/2,N);
x_h =  zeros(L,N);


% average FFT spectrum
for j = 1:N
    for i = 1:FrameN
        x_h(:,j) = xAAA((i-1)*L+1 : L*i, j);
        x_h(:,j) = x_h(:,j).*h;
        x_ps_temp = abs(fft(x_h(:,j)));
        x_ps_temp = x_ps_temp(1:L/2);
        x_ps(:,j) = x_ps(:,j) + x_ps_temp; 
    end
end



%normalize and get rid off offset of power spectrum
x_ps = 20*log10(x_ps);
offset = repmat(scaleMean, 1, N);
x_ps = x_ps - offset;



%normalized std for each sample
% for i = 1:N
%     x_ps(:,i) = x_ps(:,i)./ sqrt(sum(x_ps(:,i).*x_ps(:,i))); 
%     std = sqrt(sum(x_ps(:,i).*x_ps(:,i))); 
% end

% normalized std for each feature axis
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./ sqrt(sum(x_ps(i,:).*x_ps(i,:))); 
end

%PCA transpose
plotN = 10;
%transT = pcs'*x_ps;
transT = (V'*x_ps);


for i = 1:num_temp0
hold on
plot(transT(1:plotN ,i))
end

for i = 1:num_temp1
plot(transT(1:plotN ,i+num_temp0), 'g')
end

if(num_temp2>0)
    for i = 1:num_temp2
        plot(transT(1:plotN ,i+num_temp0+num_temp1), 'r')
    end
end

if(num_temp3>0)
    for i = 1:num_temp3
        plot(transT(1:plotN ,i+num_temp0+num_temp1+num_temp2), 'y')
    end
end

if(num_temp4>0)
    for i = 1:num_temp4
        plot(transT(1:plotN ,i+num_temp0+num_temp1+num_temp2 + +num_temp3), 'k')
    end
end
