% GS-Hardware-in-the-Loop-Ready: 平顶光生成与 MATLAB 内部实拍捕获 (闭环基石版)
clear; clc; close all;
%% --- 强制重置底层硬件 ---
fprintf('--> 正在清理后台进程，准备独占相机...\n');
imaqreset;  
pause(1.5); 

%% 0. 载入外部库文件与 DLL 初始化
addpath(genpath('SLM_Library')); 
if ~libisloaded('SecondDll')
    loadlibrary('SecondDll.dll', 'SecondDll.h');
end
temp_slm_path = fullfile(pwd, 'temp_phase_verify.bmp');

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
shift_x = 100;     % X向偏移 (避开零级光)
shift_y = 100;     % Y向偏移
target_width_um  = target_w_px * cam_pixel_pitch;
target_height_um = target_h_px * cam_pixel_pitch;
target_w_algo = round((target_width_um / 1000) / dp_focal);
target_h_algo = round((target_height_um / 1000) / dp_focal);

%% 3. 相机长连接：抓取零级光基准 (已加入曝光硬保护机制)
fprintf('--> 启动相机并保持长连接，抓取零级焦点...\n');
try
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    triggerconfig(vid, 'manual');
    
    % 获取相机支持的曝光极限，给没开光的情况加上安全锁
    try
        exp_info = propinfo(src, 'ExposureTime');
        EXP_MAX = min(exp_info.ConstraintValue(2), 1000000); 
        EXP_MIN = max(exp_info.ConstraintValue(1), 10);      
    catch
        EXP_MAX = 1000000;
        EXP_MIN = 10;
    end
    
    src.ExposureTime = 2000; 
    start(vid); pause(0.5); 
    
    % 让 SLM 黑屏以寻找纯净的零级光
    calllib('SecondDll', 'saShowImageFromFilePath', 'Black.bmp', 0, 1920, 0, 1920, 1080, 1); 
    pause(0.5);
    
    % 第一阶段自动曝光 (找焦点)
    for attempt = 1:15
        img_test = getsnapshot(vid);
        current_max = max(max(medfilt2(double(img_test), [3 3]))); 
        
        if current_max >= 250
            new_exp = round(src.ExposureTime * 0.6);
            src.ExposureTime = max(new_exp, EXP_MIN); 
            pause(0.1); 
        elseif current_max < 200
            ratio = min(240/max(current_max,1), 2.0);
            new_exp = round(src.ExposureTime * ratio);
            src.ExposureTime = min(new_exp, EXP_MAX); 
            pause(0.1); 
        else
            break; 
        end
    end
    zero_order_img = img_test;
    
    [~, max_idx] = max(zero_order_img(:));
    [cy, cx] = ind2sub(size(zero_order_img), max_idx);
    fprintf('--> 焦点坐标锁定: (X: %d, Y: %d) | 零级光曝光: %.0f us\n', cx, cy, src.ExposureTime);
    
catch ME
    error('相机独占失败！请检查是否关闭了原厂软件，或拔插相机网线/USB线重试。\n报错信息: %s', ME.message);
end

%% 4. 计算 GS 相位并打入 SLM
fprintf('--> 正在通过 GS 算法计算相位图 (50次)...\n');
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
phi_2pi = mod(phi_slm, 2*pi); 
phi_gray_1080 = uint8((phi_2pi / (2*pi)) * 255);
phi_gray_export = zeros(1080, 1920, 'uint8');
start_x = round((1920 - 1080) / 2) + 1; 
phi_gray_export(:, start_x:start_x+1079) = phi_gray_1080;
imwrite(phi_gray_export, temp_slm_path, 'bmp');

% 驱动 SLM
calllib('SecondDll', 'saShowImageFromFilePath', temp_slm_path, 0, 1920, 0, 1920, 1080, 1);
fprintf('--> 相位已打入 SLM，等待液晶偏转...\n');
pause(1.0); % 给液晶充分的物理响应时间

%% 5. 相机唤醒：MATLAB 内捕获实拍平顶光
fprintf('--> 唤醒相机，开始拉高曝光捕捉平顶光...\n');
try
    % 创建中心零级光的物理屏蔽区 (掩膜半径 150)
    [cam_H, cam_W] = size(zero_order_img); 
    [CamX, CamY] = meshgrid(1:cam_W, 1:cam_H);
    mask_zero = ((CamX - cx).^2 + (CamY - cy).^2) > 150^2; 
    
    % 起步拉高 30 倍曝光
    src.ExposureTime = min(round(src.ExposureTime * 30), EXP_MAX); 
    pause(0.2);
    
    for attempt = 1:25
        img_test = getsnapshot(vid);
        img_masked = double(img_test) .* mask_zero; % 忽略零级光测算亮度
        current_max = max(img_masked(:)); 
        
        if current_max >= 250
            src.ExposureTime = max(round(src.ExposureTime * 0.7), EXP_MIN); pause(0.1); 
        elseif current_max < 180
            ratio = min(220 / max(current_max, 1), 3.0); 
            src.ExposureTime = min(round(src.ExposureTime * ratio), EXP_MAX); pause(0.1); 
        else
            break;
        end
    end
    
    snapshot_flattop = getsnapshot(vid);
    
    % 【修改核心点】：把最终曝光数值作为普通变量存下来，防止硬件释放后找不到
    final_exposure = src.ExposureTime; 
    fprintf('--> 实拍平顶光捕获成功！最终曝光: %.0f us\n', final_exposure);
    
catch ME
    warning('平顶光实拍失败: %s', ME.message);
    snapshot_flattop = zeros(size(zero_order_img), 'uint8');
    final_exposure = 0; % 给个默认值防止画图报错
end

%% 6. 安全停机与资源释放
if exist('vid', 'var') && isvalid(vid)
    stop(vid); delete(vid); clear src vid; 
end
% 因为整个流程已经跑完，最后再 imaqreset 保证下次运行不卡死
imaqreset; 
fprintf('--> MATLAB 相机句柄已安全释放。\n');

%% 7. 终极结果对齐展示 (为后续闭环比对做准备)
figure('Name', '纯 MATLAB 闭环实验底座', 'Position', [100, 100, 1600, 500], 'Color', 'w');

% 图 1: 理论预期
U_sim = fftshift(fft2(ifftshift( Amp_G .* exp(1i * phi_slm) )));
subplot(1, 3, 1); 
imshow(abs(U_sim).^2, []); colormap(gca, 'parula');
xlim([N/2 - 200, N/2 + 200]); ylim([N/2 - 200, N/2 + 200]);
title(sprintf('1. 理论平顶光 (%dx%d)', target_w_px, target_h_px));

% 图 2: 零级光打底图
subplot(1, 3, 2); 
imshow(zero_order_img, []); colormap(gca, 'parula');
hold on; plot(cx, cy, 'r+', 'MarkerSize', 10, 'LineWidth', 2); hold off;
title('2. 实拍零级基准');

% 图 3: MATLAB 内置高曝光实拍平顶光
subplot(1, 3, 3); 
% 用高对比度 'hot' 伪彩显示
valid_pixels = double(snapshot_flattop(mask_zero)); 
display_max = 255;
if ~isempty(valid_pixels)
    sorted_vals = sort(valid_pixels(:));
    display_max = max(sorted_vals(max(round(length(sorted_vals)*0.995), 1)), 5); 
end
imshow(snapshot_flattop, [0, display_max]); colormap(gca, 'hot');
hold on; 
viscircles([cx, cy], 150, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 0.5);
text(cx, cy+180, '零级被强制屏蔽', 'Color', 'w', 'HorizontalAlignment', 'center');

% 【修改核心点】：这里调用刚才存下来的变量 final_exposure
title(sprintf('3. 实拍捕获 (曝光: %.0f us)', final_exposure));

fprintf('--> 流程结束！图像变量 snapshot_flattop 已在工作区就绪，可随时接入闭环代码！\n');