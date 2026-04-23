% GS-Hardware-in-the-Loop-Ready: 平顶光闭环反馈自适应优化版 (纯图像形态学追踪)
clear; clc; close all;
solution_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(solution_dir);
slm_library_dir = fullfile(project_root, 'SLM_Library');
if exist(slm_library_dir, 'dir')
    addpath(genpath(slm_library_dir));
else
    addpath(genpath('SLM_Library'));
end
black_bmp_path = fullfile(project_root, 'Black.bmp');
if ~exist(black_bmp_path, 'file')
    black_bmp_path = 'Black.bmp';
end
%% --- 强制重置底层硬件 ---
fprintf('--> 正在清理后台进程，准备独占相机...\n');
imaqreset;  
pause(1.5); 
%% 0. 载入外部库文件与 DLL 初始化
% SLM library path is configured above.
if ~libisloaded('SecondDll')
    loadlibrary('SecondDll.dll', 'SecondDll.h');
end
temp_slm_path = fullfile(solution_dir, 'temp_phase_verify.bmp');
%% 1. 物理参数与网格设定
N = 1080;               
p_slm = 0.0064;         
lambda = 0.000532;      
f = 300;                
D = N * p_slm;    
omega = D / 2;    
x = linspace(-D/2, D/2, N);
y = linspace(-D/2, D/2, N);
[X, Y] = meshgrid(x, y);
dp_focal = (lambda * f) / (N * p_slm); 
cam_pixel_pitch = 3.45; 
%% 2. 【核心：绝对尺寸与象限设定】
target_w_px = 200; % 平顶光宽度 (像素)
target_h_px = 150; % 平顶光高度 (像素)
shift_x = 100;     % X向偏移 (算法中向右)
shift_y = 100;     % Y向偏移 (算法中向下，即右下角)
target_width_um  = target_w_px * cam_pixel_pitch;
target_height_um = target_h_px * cam_pixel_pitch;
target_w_algo = round((target_width_um / 1000) / dp_focal);
target_h_algo = round((target_height_um / 1000) / dp_focal);
%% 3. 相机长连接：抓取零级光基准
fprintf('--> 启动相机并保持长连接，抓取零级焦点...\n');
try
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    triggerconfig(vid, 'manual');
    
    try
        exp_info = propinfo(src, 'ExposureTime');
        EXP_MAX = min(exp_info.ConstraintValue(2), 1000000); 
        EXP_MIN = max(exp_info.ConstraintValue(1), 10);      
    catch
        EXP_MAX = 1000000; EXP_MIN = 10;
    end
    
    src.ExposureTime = 2000; 
    start(vid); pause(0.5); 
    
    calllib('SecondDll', 'saShowImageFromFilePath', black_bmp_path, 0, 1920, 0, 1920, 1080, 1); 
    pause(0.5);
    
    for attempt = 1:15
        img_test = getsnapshot(vid);
        current_max = max(max(medfilt2(double(img_test), [3 3]))); 
        if current_max >= 250
            src.ExposureTime = max(round(src.ExposureTime * 0.6), EXP_MIN); pause(0.1); 
        elseif current_max < 200
            src.ExposureTime = min(round(src.ExposureTime * min(240/max(current_max,1), 2.0)), EXP_MAX); pause(0.1); 
        else
            break; 
        end
    end
    zero_order_img = img_test;
    [~, max_idx] = max(zero_order_img(:));
    [cy, cx] = ind2sub(size(zero_order_img), max_idx);
    fprintf('--> 焦点坐标锁定: (X: %d, Y: %d) | 零级光曝光: %.0f us\n', cx, cy, src.ExposureTime);
catch ME
    error('相机独占失败: %s', ME.message);
end
%% 4. 计算纯算法 GS 相位打底 (50次)
fprintf('--> 正在通过 GS 算法计算初始相位图 (纯算法迭代 50 次)...\n');
I_G = exp(-2 * (X.^2 + Y.^2) / omega^2); 
Amp_G = sqrt(I_G); 
I_T = zeros(N, N);
row_start = round(N/2 - target_h_algo/2) + shift_y + 1;
row_end   = round(N/2 + target_h_algo/2) + shift_y;
col_start = round(N/2 - target_w_algo/2) + shift_x + 1;
col_end   = round(N/2 + target_w_algo/2) + shift_x;
I_T(row_start:row_end, col_start:col_end) = 1; 
I_T = I_T * (sum(I_G(:)) / sum(I_T(:))); 
Amp_T = sqrt(I_T); 
phi_slm = 2 * pi * rand(N, N); 
for k = 1:50
    U_focal = fftshift(fft2(ifftshift( Amp_G .* exp(1i * phi_slm) )));
    phi_slm = angle(fftshift(ifft2(ifftshift( Amp_T .* exp(1i * angle(U_focal)) ))));
end
% 初始相位打入 SLM
phi_2pi = mod(phi_slm, 2*pi); 
phi_gray_export = zeros(1080, 1920, 'uint8');
start_x = round((1920 - 1080) / 2) + 1; 
phi_gray_export(:, start_x:start_x+1079) = uint8((phi_2pi / (2*pi)) * 255);
imwrite(phi_gray_export, temp_slm_path, 'bmp');
calllib('SecondDll', 'saShowImageFromFilePath', temp_slm_path, 0, 1920, 0, 1920, 1080, 1);
pause(1.0); 
%% 5. 相机唤醒：拉高曝光准备闭环
fprintf('--> 唤醒相机，自适应拉高曝光...\n');
[cam_H, cam_W] = size(zero_order_img); 
[CamX, CamY] = meshgrid(1:cam_W, 1:cam_H);
mask_zero = ((CamX - cx).^2 + (CamY - cy).^2) > 150^2; 
src.ExposureTime = min(round(src.ExposureTime * 30), EXP_MAX); pause(0.2);
for attempt = 1:25
    img_test = double(getsnapshot(vid)) .* mask_zero;
    current_max = max(img_test(:)); 
    if current_max >= 250
        src.ExposureTime = max(round(src.ExposureTime * 0.7), EXP_MIN); pause(0.1); 
    elseif current_max < 180
        src.ExposureTime = min(round(src.ExposureTime * min(220 / max(current_max, 1), 3.0)), EXP_MAX); pause(0.1); 
    else
        break;
    end
end
final_exposure = src.ExposureTime;
initial_flattop = getsnapshot(vid); % 存一张闭环前的初始状态做对比
fprintf('--> 曝光准备完毕: %.0f us。正式进入硬件闭环反馈！\n\n', final_exposure);
%% 6. 【核心】硬件在环迭代 (Camera-in-the-Loop)
fprintf('================ 开始实验实拍反馈迭代 ================\n');
closed_loop_iters = 15;   % 实验闭环次数
alpha = 0.3;              % 反馈步长 (调小至0.3防震荡)
A_weight = Amp_T;         % 初始化权重矩阵
Target_Amp_ROI = Amp_T(row_start:row_end, col_start:col_end); % 理想平顶的局部振幅

% 6.1 Use Fourier-physics mapping to define the camera ROI
% The four-point calibration file is only a validation artifact; the loop ROI
% is computed from lambda, f, SLM pitch, camera pitch, and algorithm shift.
algo_roi_h = row_end - row_start + 1;
algo_roi_w = col_end - col_start + 1;

[physics_roi, physics_info] = computePhysicsFlatTopROI( ...
    [cx, cy], ...
    [cam_H, cam_W], ...
    [shift_x, shift_y], ...
    [algo_roi_w, algo_roi_h], ...
    lambda, f, N, p_slm, cam_pixel_pitch);

cam_cx = physics_roi.center_x;
cam_cy = physics_roi.center_y;
cam_r_start = physics_roi.r_start;
cam_r_end   = physics_roi.r_end;
cam_c_start = physics_roi.c_start;
cam_c_end   = physics_roi.c_end;
cam_roi_h = physics_roi.height;
cam_roi_w = physics_roi.width;
target_mask_cam = physics_roi.mask;

fprintf('--> Physics ROI mapping locked.\n');
fprintf('    scale_ratio = %.6f camera px / algorithm px\n', physics_info.scale_ratio);
fprintf('    predicted center: (X: %d, Y: %d)\n', cam_cx, cam_cy);
fprintf('    predicted ROI: rows [%d, %d], cols [%d, %d], size = %d x %d\n', ...
        cam_r_start, cam_r_end, cam_c_start, cam_c_end, cam_roi_w, cam_roi_h);
if physics_roi.is_clipped
    warning('Physics ROI is clipped by the camera boundary. Check zero-order center or shift settings.');
end
for loop = 1:closed_loop_iters
    % 1. 实拍获取当前光场
    img_raw = double(getsnapshot(vid));
    
    % 2. 截取实拍 ROI 
    cam_crop = img_raw(cam_r_start:cam_r_end, cam_c_start:cam_c_end);
    
    % ==========================================================
    % 【光学倒像防发散机制】
    % ==========================================================
    cam_crop = rot90(cam_crop, 2);  % 根据上一轮的经验，这里必须转180度
    cam_crop = imgaussfilt(cam_crop, 1.5); % 滤除散斑
    % ==========================================================
    
    % 3. 缩放为算法尺寸
    cam_resized = imresize(cam_crop, [algo_roi_h, algo_roi_w], 'bicubic');
    
    % 4. 提取实拍振幅并归一化能量
    A_cam_roi = sqrt(max(cam_resized, 0));
    A_cam_norm = A_cam_roi * (sum(Target_Amp_ROI(:)) / sum(A_cam_roi(:)));
    
    % 计算误差 RMS (供监控)
    err_rms = sqrt(mean((Target_Amp_ROI(:) - A_cam_norm(:)).^2));
    fprintf('  Loop %02d/%02d | 实拍光斑 RMS 误差: %.4f\n', loop, closed_loop_iters, err_rms);
    
    % 5. 更新自适应权重 (Weighted GS)
    Weight_ROI = A_weight(row_start:row_end, col_start:col_end);
    Weight_ROI = Weight_ROI + alpha * (Target_Amp_ROI - A_cam_norm);
    Weight_ROI(Weight_ROI < 0) = 0; % 防止过度补偿出现负值
    A_weight(row_start:row_end, col_start:col_end) = Weight_ROI;
    
    % 6. 利用新权重进行少量的内部数字 GS 迭代
    for k = 1:5
        U_focal = fftshift(fft2(ifftshift( Amp_G .* exp(1i * phi_slm) )));
        phi_slm = angle(fftshift(ifft2(ifftshift( A_weight .* exp(1i * angle(U_focal)) ))));
    end
    
    % 7. 下发新相位到 SLM
    phi_2pi = mod(phi_slm, 2*pi); 
    phi_gray_export(:, start_x:start_x+1079) = uint8((phi_2pi / (2*pi)) * 255);
    imwrite(phi_gray_export, temp_slm_path, 'bmp');
    calllib('SecondDll', 'saShowImageFromFilePath', temp_slm_path, 0, 1920, 0, 1920, 1080, 1);
    
    pause(0.6); % 等待液晶偏转完毕，准备下一次拍照
end
% 获取闭环优化后的最终结果
final_flattop = getsnapshot(vid);
fprintf('================ 闭环优化结束 ================\n');
%% 7. 安全停机与资源释放
if exist('vid', 'var') && isvalid(vid)
    stop(vid); delete(vid); clear src vid; 
end
imaqreset; 
fprintf('--> MATLAB 相机句柄已安全释放。\n');
%% 8. 闭环前后对比展示
figure('Name', '纯算法 VS 硬件闭环优化对比', 'Position', [100, 100, 1600, 500], 'Color', 'w');
% 图 1: 算法理想预测
U_sim = fftshift(fft2(ifftshift( Amp_G .* exp(1i * phi_slm) )));
subplot(1, 3, 1); 
imshow(abs(U_sim).^2, []); colormap(gca, 'parula');
xlim([N/2 - 200, N/2 + 200]); ylim([N/2 - 200, N/2 + 200]);
title('1. 算法层面: 理论光斑');
% 图 2: 闭环前的实拍 (仅纯算法打底)
subplot(1, 3, 2); 
imshow(initial_flattop, [0, 255]); colormap(gca, 'hot');
hold on; 
viscircles([cx, cy], 150, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 0.5);
% 画出自动寻找的 ROI 抓取框
rectangle('Position', [cam_c_start, cam_r_start, cam_roi_w, cam_roi_h], 'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '-.');
title('2. 实拍: 闭环前 (仅GS算法)');
% 图 3: 闭环后的实拍
subplot(1, 3, 3); 
imshow(final_flattop, [0, 255]); colormap(gca, 'hot');
hold on; 
viscircles([cx, cy], 150, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 0.5);
rectangle('Position', [cam_c_start, cam_r_start, cam_roi_w, cam_roi_h], 'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '-.');
title('3. 实拍: 闭环反馈 15 次后');