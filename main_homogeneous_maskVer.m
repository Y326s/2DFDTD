clear;
n = 1200;
epsilon_0 = 8.85*(10^(-12));
miu_0 = 4*pi()*10^(-7);

Chara_f = 298.428*(10^6); % Characteristic frequency, Hz
Chara_T = 1/Chara_f;
Max_f = 3.351*Chara_f;

% n for t grides, i,j for x,y grids
% Cell center size: 300 * 500 * n
% For Ez: 301 * 501 * n
% For Hx: 300 * 501 * n
% For Hy: 301 * 500 * n
% CAUTION!! above size in (j, i, n) order 
Ez = zeros(301, 501, n);
Hx = zeros(300, 501, n);
Hy = zeros(301, 500, n);
Jz_filter = zeros(299, 499, n);

cell_mask = zeros(300, 500);
cell_mask(1:100,1:300) = 1;
cell_mask(100:200,50:400) = 1;
cell_mask(200:300,250:500) = 1;
imagesc(cell_mask);
colormap jet
colorbar

core_Ez = ones(2,2);
core_Hx = ones(1,2);
core_Hy = ones(2,1);
Ez_mask = imdilate(cell_mask, core_Ez, 'full');
Ez_boundary_mask = bwperim(Ez_mask);
Ez_cen_mask = Ez_mask - Ez_boundary_mask;
Hx_mask = imdilate(cell_mask, core_Hx, 'full');
Hy_mask = imdilate(cell_mask, core_Hy, 'full');



epsilon_r_matrix = ones(300, 500);
miu_r_matrix = ones(300, 500);
epsilon_matrix = epsilon_r_matrix .* epsilon_0;
miu_matrix = miu_r_matrix .* miu_0;

c_matrix = 1./sqrt(epsilon_matrix.*miu_matrix);
c_min = min(c_matrix,[],'all');
c_max = max(c_matrix,[],'all');
lambda_min = c_min/Max_f; %

delta_x = lambda_min/12;
delta_y = lambda_min/12;
delta_t = lambda_min/(24*c_max); % To be confirmed 

%
%epsilon_matrix = epsilon_r_matrix .* epsilon_0;
%miu_matrix = miu_r_matrix .* miu_0;

core_miu_x = [0.5, 0.5];
core_miu_y = [0.5; 0.5];
core_epsilon = [0.25,0.25;0.25,0.25];
miu_for_Hx = conv2(miu_matrix, core_miu_x);
miu_for_Hy = conv2(miu_matrix, core_miu_y);
epsilon_for_Ez = conv2(epsilon_matrix, core_epsilon);

miu_matrix_Hx = [miu_matrix(:,1), miu_for_Hx(:,2:end-1), miu_matrix(:,end)];
miu_matrix_Hy = [miu_matrix(1,:); miu_for_Hy(2:end-1,:); miu_matrix(end,:)];
epsilon_matrix = epsilon_for_Ez(2:end-1, 2:end-1);

Ez_x_diff = zeros(301, 500, n);
Ez_y_diff = zeros(300, 501, n);
Hy_x_diff = zeros(301, 499, n);
Hx_y_diff = zeros(299, 501, n);



source_t_range = floor((Chara_T-delta_t)/delta_t) + 1

for ss = 2:source_t_range
    Jz_filter(49,49,ss) = Jz_BHW((ss-3/2)*delta_t,Chara_T); % To be confirmed
    %49,49
end
Jz_filter(49,49,150)



for pp=2:n
    Hx(:,:,pp) = Hx(:,:,pp-1) + (delta_t./miu_matrix_Hx) .* (Ez_y_diff(:,:,pp-1)) ./ delta_y;
    Hx = Hx .* Hx_mask;
    Hy(:,:,pp) = Hy(:,:,pp-1) + (delta_t./miu_matrix_Hy) .* (Ez_x_diff(:,:,pp-1)) ./ delta_x;
    Hy = Hy .* Hy_mask;
    Hy_x_diff(:,:,pp) = diff(Hy(:,:,pp),1,2);
    Hx_y_diff(:,:,pp) = diff(Hx(:,:,pp));
    
    Ez(2:end-1,2:end-1,pp) = Ez(2:end-1,2:end-1,pp-1) + ...
        (delta_t./epsilon_matrix) .* (Hy_x_diff(2:end-1,:,pp)./delta_x + ...
        Hx_y_diff(:,2:end-1,pp)./delta_y) + (delta_t./epsilon_matrix) .* Jz_filter(:,:,pp); %%% + !!!
    Ez = Ez .* Ez_cen_mask;
    Ez_y_diff(:,:,pp) = diff(Ez(:,:,pp));
    Ez_x_diff(:,:,pp) = diff(Ez(:,:,pp),1,2);
end

% save('./EHresult', 'Ez', 'Hx', 'Hy', 'n')


















