clc; clear all; close all;
arr = load('data_d3.mat', '-ascii');
arr = arr;
clut_rmv = arr(:,[1:107]);
clut_rmv(:,[1:20]) = 0;
clut_rmv(:,[21:24]) = 0.8*clut_rmv(:,[21:24]);
% range=0:62.5/256:62.5-62.5/256;
range1=linspace(0,62.5,256);
range = range1([1:107]);
% range = 0:25/100:26;
dopp = 1:length(arr(:,1));
dopp_g = dopp*180/length(dopp);
figure(2); surfl(range,dopp_g,abs(clut_rmv));
title('3D plot of rotating AWR1843 - Human at 5m-60deg, 7m-90deg and 9m-120deg');
grid on; grid minor; xlabel('Range (m)'); ylabel('angle (degree)'); zlabel('Amplitude')

% figure(1); imagesc(range,dopp_g,mag2db(mean_signal_signal_1_cut')); colorbar;
figure(3); imagesc(range,dopp_g,mag2db(abs(clut_rmv))); colorbar;
% title('3D plot of rotating AWR1843 - Case (a)');
grid on; grid minor; xlabel('Range (m)'); ylabel('angle (degree)'); zlabel('Amplitude');

%% CFAR implementation
Tr = 15; Td = 15; % Training cells
% Guard Cells in both dimensions around the Cell under test (CUT)
Gr = 10; Gd = 10;
% offset the threshold by SNR value in dB
offset = 5.5;
TGr = Tr + Gr;
TGd = Td + Gd;
RDM = clut_rmv;
result = zeros(size(RDM,1), size(RDM,2));
factor = 10;
for i = TGr + 1 : size(RDM,1) - (TGr)
    for j = TGd + 1 : size(RDM,2) - (TGd)
        
        A = RDM((i-TGr):(i+TGr), (j-TGd):(j+TGd));
%         A = db2pow(A);
        
        A((Tr+1):(Tr + 2*Gr +1), (Td+1):(Td + 2*Gd +1)) = 0.0;
        
        Asum = sum(A,'all');
        noise_level = Asum / (size(A,1)*size(A,2) - (2*Gr+1)*(2*Gd+1));
        noise_level = pow2db(noise_level) + factor;
        
        if RDM(i,j) > noise_level
            result(i,j) = RDM(i,j);
        else
            result(i,j) = 0;
        end
    end
end

figure('Name','CA-CFAR Filtered RDM'); 
surf(range,dopp_g,mag2db(result)); colorbar;
xlabel('Freq (Hz)'), ylabel('Range (m)'); axis tight; grid on; grid minor;



pos_5 = clut_rmv(:,21); [max_5 poss_5] = max(pos_5); angle_5 = dopp_g(poss_5);
pos_7 = clut_rmv(:,29); [max_7 poss_7] = max(pos_7); angle_7 = dopp_g(poss_7); %30
pos_9 = clut_rmv(:,38); [max_9 poss_9] = max(pos_9); angle_9 = dopp_g(poss_9);
pos_11 = clut_rmv(:,47);[max_11 poss_11] = max(pos_11); angle_11 = dopp_g(poss_11);
pos_13 = clut_rmv(:,55);[max_13 poss_13] = max(pos_13); angle_13 = dopp_g(poss_13);%
pos_15 = clut_rmv(:,63);[max_15 poss_15] = max(pos_15); angle_15 = dopp_g(poss_15);
pos_17 = clut_rmv(:,71);[max_17 poss_17] = max(pos_17); angle_17 = dopp_g(poss_17);
pos_19 = clut_rmv(:,79);[max_19 poss_19] = max(pos_19); angle_19 = dopp_g(poss_19);
pos_21 = clut_rmv(:,88);[max_21 poss_21] = max(pos_21); angle_21 = dopp_g(poss_21);%87
pos_23 = clut_rmv(:,95);[max_23 poss_23] = max(pos_23); angle_23 = dopp_g(poss_23);
pos_25 = clut_rmv(:,104);[max_25 poss_25] = max(pos_25); angle_25 = dopp_g(poss_25);%103



mean_pos_9 = mean(pos_9); std_pos_9 = std(pos_9); var_pos_9 = var(pos_9);
mean_pos_7 = mean(pos_7); std_pos_7 = std(pos_7); var_pos_7 = var(pos_7);
mean_pos_5 = mean(pos_5); std_pos_5 = std(pos_5); var_pos_5 = var(pos_5);
area_5 = abs(trapz(pos_5-mean(pos_5))); area_7 = abs(trapz(pos_7-mean(pos_7))); 
area_9 = abs(trapz(pos_9-mean(pos_9)));
figure(2); plot(dopp_g,pos_5); hold on; plot(dopp_g,pos_7); hold on; plot(dopp_g,pos_9); 
legend('pos 5','pos 7','pos 9'); xlabel('Angle (degree)'); ylabel('Amplitude');
grid on; grid minor;
title(['mean_5=' num2str(mean_pos_5) ', std_5=' num2str(std_pos_5) ',var_5='...
    num2str(var_pos_5) ', area_5=' num2str(area_5) ', mean_7=' num2str(mean_pos_7) ', std_7=' num2str(std_pos_7) ',var_7='...
    num2str(var_pos_7) ', area_7=' num2str(area_7) ', mean_9=' num2str(mean_pos_9) ', std_9=' num2str(std_pos_9) ',var_9='...
    num2str(var_pos_9) ', area_9=' num2str(area_9)],'Interpreter','tex');
% str = {['area_pkr= ' num2str(var_pos_5)],['var_pkr= ' num2str(var_pos_5)]};
% str = strcat('var_pkr=',num2str(var_pos_5), ' area_pkr=',num2str(var_pos_5), ' height_pkr=',num2str(var_pos_5));
% annotation('textbox', [0.179 0.849 0.2 0.037], 'string', str,'FitBoxToText','on','Fontsize',12,'Interpreter','tex');
% annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');




% Define sample data consisting of floating point numbers
m = clut_rmv;
% Sort the array
sortedm = sort(m(:), 'descend');
% Find the 10% index
index10 = ceil(0.03 * numel(m));
% Determine the value that 10% of the values lie above
value10 = sortedm(index10);
% Get binary "map" of where these highest 10% of values live:
binaryImage = m >= value10;
% Extract a matrix with only the highest 10% of floating point values:
output = m .* binaryImage;
figure(3);
surf(range,dopp_g,output); colorbar;
v5=[]; p5=[];angle_new=[]; range_new=[]; p5_new=[];v5_new=[];
for i = 1:length(output(:,1))
    [v p] = max(output(i,:)); 
    v5=[v5 v]; p5=[p5 p];
    [v5_n p5_n] = max(output(:,p));
    v5_new=[v5_new v5_n]; p5_new=[p5_new p5_n];
    angle5_new = dopp_g(p5_n); range5_new = range(p);
    angle_new = [angle_new angle5_new]; range_new=[range_new range5_new];
end

p=0; u=1; ang=0;
for t=2:output(1,:)
    diff = abs(angle_new(t) - angle_new(t-1));
    if diff<=5
        p = p + 1;
        ang(u) = angle_new(t);
    else
        p = p;
        u = u+1;
    end
end
v = nonzeros(ang');
final_angle= unique(v);    
hh=[]; final_range=[];

for i=1:length(final_angle)
    hh(i) = find(dopp_g == final_angle(i));
    final_range(i) = range_new(hh(i));
end
