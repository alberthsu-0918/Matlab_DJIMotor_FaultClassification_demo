% clc;
% MCM Test for compressor.  
% After run 'pcaTest_MCM_fDomain.m' to get pcs, the transpose matrix
% Run this m file to see if the pcs can project other raw data to a better
% 3d axis.
% Further clustering is needed
% 2017/11/23
%clear all;

L = 4096; %load 4096 samples of wavefile
%N = 20; %number of oberservation
N = 32;
FrameN = 5;
num_temp0 = 0;
num_temp1 = 0;
num_temp2 = 0;
num_temp3 = 0;
num_temp4 = 0;



%Dir0 = 'C:\Matlab\work\PCA\MCM\2017_11_20\';
Dir0 = 'C:\Matlab\work\PCA\MCM\sound2017_08_15_all\';
temp0=dir([Dir0,'*.csv']);
num_temp0=length(temp0);
%
% %Dir1 = 'C:\Matlab\work\PCA\MCM\2017_10_03\'; %2017_0603_bfFail\ %2017_06_04_TheBreakdingDAy\
% Dir1 = 'C:\Matlab\work\PCA\MCM\sound2017_10_17\';
% temp1=dir([Dir1,'*.csv']);
% num_temp1=length(temp1);
% % % % 
% Dir2 = 'C:\Matlab\work\PCA\MCM\sound2017_09_05\'; %2017_08_14to19
% temp2=dir([Dir2,'*.csv']);
% num_temp2=length(temp2);
% % % % 
% Dir3 = 'C:\Matlab\work\PCA\MCM\sound2017_08_15\'; %2017_07_03
% temp3=dir([Dir3,'*.csv']);
% num_temp3=length(temp3);
% %
% Dir4 = 'C:\Matlab\work\PCA\MCM\sound2017_06_20\'; %2017_0603_bfFail\ %2017_06_04_TheBreakdingDAy\
% temp4=dir([Dir4,'*.csv']);
% num_temp4=length(temp4);

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

%x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32];



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
% using pre-built pca transfer model's axis center
% Therefore new data need to minus pre-built offset, scaleMean.



x_ps = 20*log10(x_ps);
offset = repmat(scaleMean, 1, N);
x_ps = x_ps - offset;

% normalized std for each sample
% for i = 1:N
%     x_ps(:,i) = x_ps(:,i)./ sqrt(sum(x_ps(:,i).*x_ps(:,i))); 
% end

% normalized std for each spectrum data
% for i = 1:L/2
%     x_ps(i,:) = x_ps(i,:)./ sqrt(sum(x_ps(i,:).*x_ps(i,:))); 
% end


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



