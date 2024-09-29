% clear;

% n for t grides, i,j for x,y grids
%
% Ez size: 300 * 500 * n
% Hx size: 299 * 500 * n
% Hy size: 300 * 499 * n
% CAUTION!! above size in (j, i, n) order
%
%
%
% convolution！！
%
%
%

%--------------------------CONSTANTS----------------------------------
epsilon_0 = 8.85*(10^(-12));
miu_0 = 4*pi()*10^(-7);

Chara_f = 298.428*(10^6); % Characteristic frequency, Hz
Chara_T = 1/Chara_f;
Max_f = 3.351*Chara_f;

%--------------------------INITIATE-----------------------------------

mask_cell = zeros(300, 500);
mask_cell(1:100,1:300) = 1;
mask_cell(100:200,50:400) = 1;
mask_cell(200:300,250:500) = 1;
imagesc(mask_cell);
colormap jet
colorbar










n = 1000;
Ez = zeros(300, 500, n);
Hx = zeros(299, 500, n);
Hy = zeros(300, 499, n);

Jz_filter = zeros(298, 498, n);

epsilon_r_matrix = ones(299, 499);
miu_r_matrix = ones(299, 499);



%----------------------------MEDIUM-----------------------------------
% Set different medium inside
epsilon_r_matrix(100:200, 100:150) = 1;
miu_r_matrix(100:200, 100:150) = 1;

epsilon_matrix = epsilon_r_matrix .* epsilon_0;
miu_matrix = miu_r_matrix .* miu_0;
c_matrix = 1./sqrt(epsilon_matrix.*miu_matrix);


% Test c_min
c_min = min(c_matrix,[],'all');
% c_min = 1/sqrt(9*epsilon_0*miu_0);


c_max = max(c_matrix,[],'all');
lambda_min = c_min/Max_f; %

delta_x = lambda_min/12;
delta_y = lambda_min/12;
delta_t = lambda_min/(24*c_max); % To be confirmed 

c_min
c_max
delta_x
delta_y
delta_t
%---------------------------------------------------------------------



epsilon_x_extend = [epsilon_matrix(:,1), epsilon_matrix] + [epsilon_matrix, epsilon_matrix(:,end)]./2;
epsilon_y_extend = [epsilon_matrix(1,:); epsilon_matrix] + [epsilon_matrix; epsilon_matrix(end,:)]./2;
miu_x_extend = [miu_matrix(:,1), miu_matrix] + [miu_matrix, miu_matrix(:,end)]./2;
miu_y_extend = [miu_matrix(1,:); miu_matrix] + [miu_matrix; miu_matrix(end,:)]./2;
conv_core = [0.25,0.25;0.25,0.25];
epsilon_conv2 = conv2(conv_core, epsilon_matrix);
epsilon_conv2_center = epsilon_conv2(2:end-1,2:end-1);

% Ez_x_diff = diff(Ez,1,2);
% Ez_y_diff = diff(Ez);
% Hy_x_diff = diff(Hy,1,2);
% Hx_y_diff = diff(Hx);

Ez_x_diff = zeros(300, 499, n);
Ez_y_diff = zeros(299, 500, n);
Hy_x_diff = zeros(300, 498, n);
Hx_y_diff = zeros(298, 500, n);

% Ez = zeros(300, 500, n);
% Hx = zeros(299, 500, n);
% Hy = zeros(300, 499, n);

source_t_range = floor((Chara_T-delta_t)/delta_t) + 1



for pp = 2:source_t_range
    Jz_filter(49,49,pp) = Jz_BHW((pp-3/2)*delta_t,Chara_T); % To be confirmed
    %49,49
end
Jz_filter(49,49,150)

% for pp = 2:2
%     Jz_filter(49,49,pp) = Jz_BHW((100)*delta_t,Chara_T); % To be confirmed
%     %49,49
% end
% %Jz_filter(49,49,1:source_t_range)




for kk = 2:n
    Hx(:,:,kk) = Hx(:,:,kk-1) + (delta_t./miu_x_extend) .* (Ez_y_diff(:,:,kk-1)) ./ delta_y; %%%%% + !!!
    Hy(:,:,kk) = Hy(:,:,kk-1) + (delta_t./miu_y_extend) .* (Ez_x_diff(:,:,kk-1)) ./ delta_x;
    Hy_x_diff(:,:,kk) = diff(Hy(:,:,kk),1,2);
    Hx_y_diff(:,:,kk) = diff(Hx(:,:,kk));
    
    
    Ez(2:end-1,2:end-1,kk) = Ez(2:end-1,2:end-1,kk-1) + (delta_t./epsilon_conv2_center) .* (Hy_x_diff(2:end-1,:,kk)./delta_x + Hx_y_diff(:,2:end-1,kk)./delta_y) + (delta_t./epsilon_conv2_center) .* Jz_filter(:,:,kk); %%% + !!!
    
    Ez_y_diff(:,:,kk) = diff(Ez(:,:,kk));
    Ez_x_diff(:,:,kk) = diff(Ez(:,:,kk),1,2);

    % E time diff use same time grides
end

% test1 = Hy(40:60,40:60,2);
% test2 = Hy(40:60,40:60,3);
% test3 = Hy(40:60,40:60,4);
test4 = Ez(40:60,40:60,2);
test5 = Ez(40:60,40:60,3);
test6 = Ez(40:60,40:60,4);
% test7 = Hx(40:60,40:60,2);
% test8 = Hx(40:60,40:60,3);
% test9 = Hx(40:60,40:60,4);




% 
% Ez_2D = imagesc(Ez(:,:,500));
% colormap jet
% colorbar

for qq = 1:10:n
%     figure;
    
    
    imagesc(Ez(:,:,qq));
    xlabel(qq);
    colormap jet
    colorbar
    caxis([-0.1,0.1]);
    pause(0.3);
%     close all
%     
%     
% %     plot([1:size(Ey,1)],Ey(:,qq),'r');
% %     xlabel(['t=' num2str(delta_t*(qq-1)*1e9) 'ns']);
% %     axis([0 size(Ey,1) 1.1*min(Ey,[],'ALL') 1.1*max(Ey,[],'ALL')]);
end










