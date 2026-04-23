% =========================================================================
% 光学闭环基石：带零级光屏蔽与【物理衍射定律校验】的绝对像素对准程序
% =========================================================================
clear; clc; close all;

%% 0. 载入外部库文件与 DLL 初始化
fprintf('--> 正在初始化硬件用于光学对准...\n');
imaqreset; pause(1.5);
addpath(genpath('SLM_Library')); 
if ~libisloaded('SecondDll')
    loadlibrary('SecondDll.dll', 'SecondDll.h');
end
temp_slm_path = fullfile(pwd, 'temp_calibration.bmp');

%% 0.5 【核心配置】物理光学系统参数 (请根据实际光路修改)
% -------------------------------------------------------------------------
PHYSICS.lambda = 532e-6;      % 激光波长 (单位：毫米，例如 532 nm = 532e-6 mm)
PHYSICS.f_lens = 300;         % 傅里叶透镜/成像透镜焦距 (单位：毫米)
PHYSICS.p_slm  = 6.4e-3;      % SLM 物理像素尺寸 (单位：毫米，例如 8 um = 8e-3 mm)
PHYSICS.p_cam  = 3.45e-3;     % 相机 物理像素尺寸 (单位：毫米，例如 3.45 um = 3.45e-3 mm)
PHYSICS.N_algo = 1080;        % 算法生成网格的分辨率尺寸
% -------------------------------------------------------------------------

%% 1. 相机长连接与零级焦点定位
try
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    triggerconfig(vid, 'manual');
    src.ExposureTime = 2000; 
    start(vid);
    pause(0.5); 
    
    calllib('SecondDll', 'saShowImageFromFilePath', 'Black.bmp', 0, 1920, 0, 1920, 1080, 1);
    pause(0.5);
    
    fprintf('--> 正在自动调节相机曝光抓取零级焦点...\n');
    max_attempts = 20;       
    target_max_val = 240;    
    min_acceptable = 200;    
    
    for attempt = 1:max_attempts
        img_test = getsnapshot(vid);
        img_test_filtered = medfilt2(double(img_test), [3 3]); 
        current_max = max(img_test_filtered(:)); 
        
        if current_max >= 250
            src.ExposureTime = max(round(src.ExposureTime * 0.6), 10);
            pause(0.1); 
        elseif current_max < min_acceptable
            ratio = target_max_val / max(current_max, 1);
            ratio = min(ratio, 2.0); 
            src.ExposureTime = min(round(src.ExposureTime * ratio), 1000000);
            pause(0.1);
        else
            fprintf('--> 零级光曝光调节成功！基准曝光: %.0f us\n', src.ExposureTime);
            zero_order_img = img_test; 
            img_filtered = img_test_filtered; 
            break;
        end
        if attempt == max_attempts
            warning('曝光未收敛，采用最后一次抓图。');
            zero_order_img = img_test; img_filtered = img_test_filtered;
        end
    end
    zero_exposure = src.ExposureTime; 
    
catch ME
    error('相机独占失败: %s', ME.message);
end

%% 2. 锁定零级坐标，并生成“零级黑洞掩膜”
[max_val, max_idx] = max(img_filtered(:));
[cy, cx] = ind2sub(size(img_filtered), max_idx);
fprintf('--> 【零级光中心已锁定】: (X: %d, Y: %d)\n', cx, cy);
[cam_H, cam_W] = size(zero_order_img);
[CamX, CamY] = meshgrid(1:cam_W, 1:cam_H);
mask_zero = ((CamX - cx).^2 + (CamY - cy).^2) > 80^2; 

%% 3. 定义 4 个探测点 (算法矩阵坐标，以算法中心为基准)
N = PHYSICS.N_algo; 
[X_grid, Y_grid] = meshgrid(1:N, 1:N);
probe_dist = 80; 
algo_pts = [ probe_dist,  probe_dist; 
            -probe_dist,  probe_dist; 
            -probe_dist, -probe_dist; 
             probe_dist, -probe_dist];
cam_pts = zeros(4, 2); 
probe_colors = ['r', 'g', 'b', 'm']; 

%% 4. 发射探针并记录实拍坐标 (带独立动态曝光调节)
fprintf('--> 开始发射 4 点探针进行绝对坐标映射...\n');
figure('Name', '实拍探针捕获过程', 'Position', [100, 200, 1000, 800]);

for i = 1:4
    u = algo_pts(i, 1); 
    v = algo_pts(i, 2);
    phase_map = mod(2 * pi * (u .* X_grid / N + v .* Y_grid / N), 2 * pi);
    
    phi_gray = uint8((phase_map / (2*pi)) * 255);
    phi_export = zeros(1080, 1920, 'uint8');
    start_x = round((1920 - 1080) / 2) + 1; 
    phi_export(:, start_x:start_x+1079) = phi_gray;
    imwrite(phi_export, temp_slm_path, 'bmp');
    calllib('SecondDll', 'saShowImageFromFilePath', temp_slm_path, 0, 1920, 0, 1920, 1080, 1);
    
    pause(0.8); 
    
    fprintf('--> 正在为 探针 %d 自动搜索完美曝光...\n', i);
    src.ExposureTime = min(zero_exposure * 50, 500000); 
    
    max_probe_attempts = 15;
    target_probe_val = 230; 
    min_probe_acceptable = 200; 
    
    for att = 1:max_probe_attempts
        pause(max(0.1, src.ExposureTime / 1000000 + 0.1)); 
        
        img_temp_raw = double(getsnapshot(vid));
        img_temp_masked = img_temp_raw .* mask_zero; 
        
        img_temp_filtered = medfilt2(img_temp_masked, [3 3]);
        cur_max = max(img_temp_filtered(:)); 
        
        if cur_max >= 250
            src.ExposureTime = max(round(src.ExposureTime * 0.7), 10);
        elseif cur_max < min_probe_acceptable
            ratio = target_probe_val / max(cur_max, 1);
            ratio = min(ratio, 2.5); 
            src.ExposureTime = min(round(src.ExposureTime * ratio), 3000000); 
        else
            fprintf('    ✅ 探针 %d 曝光调节达标！峰值亮度: %d, 最终曝光: %.0f us\n', i, cur_max, src.ExposureTime);
            img_probe_raw = img_temp_raw;
            img_probe_masked = img_temp_masked;
            break;
        end
        if att == max_probe_attempts
            img_probe_raw = img_temp_raw; img_probe_masked = img_temp_masked;
        end
    end
    
    max_probe_val = max(img_probe_masked(:));
    bg_thresh = max_probe_val * 0.4; 
    img_clean = img_probe_masked; 
    img_clean(img_clean < bg_thresh) = 0; 
    
    mass = sum(img_clean(:));
    if mass > 0
        cam_cx = sum(img_clean(:) .* CamX(:)) / mass;
        cam_cy = sum(img_clean(:) .* CamY(:)) / mass;
    else
        cam_cx = cx; cam_cy = cy;
    end
    cam_pts(i, :) = [cam_cx, cam_cy];
    
    subplot(2, 2, i);
    imshow(img_probe_raw, [0 255]); colormap(gca, hot); hold on;
    viscircles([cx, cy], 80, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 0.5);
    plot(cam_cx, cam_cy, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
    text(cam_cx + 20, cam_cy, sprintf('P%d', i), 'Color', probe_colors(i), 'FontSize', 14, 'FontWeight', 'bold');
    drawnow; hold off;
end

%% 5. 恢复黑屏并释放硬件
calllib('SecondDll', 'saShowImageFromFilePath', 'Black.bmp', 0, 1920, 0, 1920, 1080, 1);
stop(vid); delete(vid); clear src vid; imaqreset;

%% 6. 几何变换映射推算
cam_pts_centered = cam_pts - [cx, cy]; 
tform = fitgeotrans(algo_pts, cam_pts_centered, 'affine');
T = tform.T; 
angle_deg = rad2deg(atan2(-T(1,2), T(1,1)));

algo_width = abs(algo_pts(1,1) - algo_pts(2,1)); 
cam_width  = sqrt((cam_pts(1,1)-cam_pts(2,1))^2 + (cam_pts(1,2)-cam_pts(2,2))^2);
actual_scale = cam_width / algo_width;

%% 7. 波光学理论验证分析
% 依据公式： u_cam = (f * lambda) / (N * P_slm * P_cam) * u_algo
theory_scale = (PHYSICS.f_lens * PHYSICS.lambda) / ...
               (PHYSICS.N_algo * PHYSICS.p_slm * PHYSICS.p_cam);

% 如果有倍率放大（比如 4f 系统有 m 倍），等效焦距为：
actual_equivalent_focal = (actual_scale * PHYSICS.N_algo * PHYSICS.p_slm * PHYSICS.p_cam) / PHYSICS.lambda;
scale_error = (actual_scale - theory_scale) / theory_scale * 100;

%% 8. 打印终极物理报告
fprintf('\n================== 像素匹配与物理标定报告 ==================\n');
fprintf('【光路几何诊断】\n');
fprintf('1. 绝对零级中心: X = %d, Y = %d\n', cx, cy);
fprintf('2. 光路相对旋转角: %.2f 度 (SLM面板到相机芯片的物理偏置)\n', angle_deg);
fprintf('\n【物理衍射诊断】\n');
fprintf('3. 理论物理放大率: 1 算法像素 = %.3f 相机像素\n', theory_scale);
fprintf('4. 实测物理放大率: 1 算法像素 = %.3f 相机像素\n', actual_scale);
fprintf('5. 尺度误差偏差率: %+.2f %%\n', scale_error);
fprintf('\n【误差来源分析建议】:\n');
if abs(scale_error) < 3
    fprintf('   - ✅ 误差极小！系统符合标准理想傅里叶透镜模型。\n');
else
    fprintf('   - ❌ 存在显著倍率差异。原因可能是：\n');
    fprintf('     a) 透镜的实际焦距与标称值(%.1f mm)有偏差，推算的真实等效焦距为 %.1f mm\n', PHYSICS.f_lens, actual_equivalent_focal);
    fprintf('     b) 如果你用了两个透镜(比如4f系统)，它们的焦距并不完全相等，引入了光学放大率 M。\n');
    fprintf('     c) 相机传感器没有严格放置在透镜的焦平面上。\n');
end
fprintf('======================================================\n');

%% 9. 数据持久化保存
save_filename = 'Calibration_Mapping_Data.mat';
calib_data.algo_pts = algo_pts;
calib_data.cam_pts_centered = cam_pts_centered;
calib_data.affine_matrix = T;
calib_data.rotation_angle_deg = angle_deg;
calib_data.actual_scale = actual_scale;
calib_data.theory_scale = theory_scale;
save(save_filename, 'calib_data');