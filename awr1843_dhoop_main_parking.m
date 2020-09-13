%%% PROGRAM TO ANALYSE DATA RECORDED USING DCA1000 AND AWR1243 for the file saved in ppro/data1.
clc; clear all; close all;

%% global variables
% Based on sensor configuration.
   numADCBits = 16; % number of ADC bits per sample.
   numADCSamples = 256; % number of ADC samples per chirp.
   numRx = 4; % number of receivers in AWR1243.
   chirpSize = numADCSamples*numRx;
   chirploops= 128; % No. of of chirp loops.  
   numLanes = 2; % do not change. number of lanes is always 4
   isReal = 0; % set to 1 if real only data, 0 if complex data.
   numFrames = 200; 
   numChirps = 1;
   sampleRate = 10; % [Msps]
   timeStep = 1/sampleRate;    % [us]
   chirpPeriod = numADCSamples * timeStep ; % [us]
   plotEnd = numADCSamples * numChirps*numFrames; %for considering all frames.
   Dx = numADCSamples * numChirps ;
   timeEnd = (plotEnd-1) * timeStep;

%% read file
% read .bin file
fid = fopen('adc_data.bin','r');
% adcData = fread(fid, 'int16');
% fclose(fid);
% fileSize = size(adcData, 1);


% % % % % % % % 

adcData = fread(fid, 'int16');
% if 12 or 14 bits ADC per sample compensate for sign extension
if numADCBits ~= 16
    l_max = 2^(numADCBits-1)-1;
    adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
end
fclose(fid);
fileSize = size(adcData, 1);
% real data reshape, filesize = numADCSamples*numChirps
if isReal
    numChirps = fileSize/numADCSamples/numRx;
    LVDS = zeros(1, fileSize);
    %create column for each chirp
    LVDS = reshape(adcData, numADCSamples*numRx, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
else
    % for complex data
    % filesize = 2 * numADCSamples*numChirps
    numChirps = fileSize/2/numADCSamples/numRx;
    LVDS = zeros(1, fileSize/2);
    %combine real and imaginary part into complex data
    %read in file: 2I is followed by 2Q
    counter = 1;
    for i=1:4:fileSize-1
        LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
        LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
        counter = counter + 2;
    end
    % create column for each chirp
    LVDS = reshape(LVDS, numADCSamples*numRx, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
end

%organize data per RX
adcData = zeros(numRx,numChirps*numADCSamples);
adc_r1=[]; adc_r2=[]; adc_r3=[]; adc_r4=[];

for i=1:256
    adc_r1(:,i) = LVDS(:,i);
    adc_r2(:,i) = LVDS(:,i+256);
    adc_r3(:,i) = LVDS(:,i+512);
    adc_r4(:,i) = LVDS(:,i+768);
end
% filename = '0_case_0.xlsx';
% xlswrite(filename,LVDS);
% fft_zero_all = xlsread('fft_zero_all.xlsx');
max_peak=3;
new_pos=[]; new_val=[]; easy_positions=zeros(400,max_peak);
% easy_positions=zeros(200,200);

for n = 1 : 400
    u = n - 1;
    ff_1_tp = LVDS(1+128*u:128+128*u,:);
            
    D = 2; N = 25;
    ff_1_tp_new = ff_1_tp';
    mean_frame = max(ff_1_tp);
    mean_frame_ch_1(i,:) = mean_frame;
    
    R_FFT = fft(ff_1_tp_new)/1024;
    R_FFT = fft(mean_frame')/1024;
    R_FFT = flip(R_FFT);
    R_FFT_new(:,n) = R_FFT;
    zero_clut_R_FFT = real(R_FFT);
%     - fft_zero_all(:,n);
    
%     zero_clut_R_FFT=fft_zero_all(:,n);
    R_FFT(1:70)=0;
    R_FFT(700:end)=0;
    R_FFT_new(:,n) = R_FFT;
    
    
    zero_clut_R_FFT(1:100) = 0;
    zero_clut_R_FFT(400:1024)=0;
%     zero_clut_R_FFT = zero_clut_R_FFT(1:512)/1024;
    zero_clut_R_FFT = zero_clut_R_FFT(1:1024)/1024;
    X2=abs(zero_clut_R_FFT);
    threshold = 1; %tolerance threshold
    X2(abs(zero_clut_R_FFT)<threshold) = 0;
    zero_clut_R_FFT=X2;
    zero_new(n,:)= zero_clut_R_FFT;
%     position calculation
    range=0:5/100:50-5/100;
%     range=0:25/512:25-25/512;
    range=0:50/1024:50-50/1024;

%         range = -25.20:0.196:24.80; 
    o=0; 
%     position calculation (easy)
    [mxpp,lopp] = max(abs(zero_clut_R_FFT));
    lopp_new(n,:)=lopp;
%     easy_positions(n) = range(lopp);
    [valpos pospos] = maxk(R_FFT, max_peak);
    pospos=pospos';
    
    easy_positions(n,:) = range(pospos);



    %Second FFT for Doppler information
    Dop = fftshift(fft(R_FFT'));
    Dop_new = Dop';
%     Dopp = [Dopp; Dop];
    [My, Ny] = size(Dop);
    Range = linspace(-N, N, Ny);
    dd = 2;
    doppler = -dd:(2*dd)/128:(dd-1/128);
%     range=0:5/100:50-5/100;
  %  subplot(2,1,1);
    plot(Range,abs(R_FFT)); 
%    plot(range,abs(zero_clut_R_FFT)); 
%     axis([0 1000 0 10^8]);
%     plot(abs(R_FFT))
%     plot(abs(Dop)); 
   drawnow;
end

dopp=1:400;
figure(3);
range=0:50/1024:50-50/1024;
surfl(dopp,range,abs(R_FFT_new));
title('3D plot of rotating- Human at 5m-60deg and 7m-90deg ')


x_circ=0;y_circ=0;r_circ=25;
    th = -pi/2:pi/50:pi/2;
    xunit = ( r_circ * cos(th) + x_circ);
    yunit = (r_circ * sin(th) + y_circ);
    figure(2);
    yyaxis left; plot(xunit, yunit);
    axis([-15 15 -25 0]); xlabel("Horizontal distance"); 
    ylabel("Vertical distance"); axis tight; grid minor; 
    title('Plot of rotating- Human at 5m-60deg and 7m-90deg')
    yyaxis right; ylim([-15 15]);
    hold on;
    % end circle
    col=['m','k','g','y','r','b','c'];
    mar=['*' 'o' '+' '-' 'd' '^' 'v' '>' '<'];
    
for i=1:400
    th1=(i-1)*pi/400;
    for j = 1:length(easy_positions(1,:))
%         if new_pos(i,j) ~= 0
        if easy_positions(i,j) ~= 0
%             r2 = abs(new_pos(i,j));
            r2 = abs(easy_positions(i,j));
%             vall = (new_val(i,j)+1)/10;
            yunit2 = r2 * cos(th1);
            xunit2 = r2 * sin(th1);
            figure(2);
            plot(xunit2,yunit2,'color', col(1),'linestyle','none','marker',mar(1),'markersize',11);
%             plot(xunit2,yunit2,'color', col(3));
            hold on;
        else
            break
        end
    end
    grid on;
    drawnow;
end