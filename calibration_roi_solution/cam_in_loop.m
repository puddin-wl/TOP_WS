% GS Hardware-in-the-Loop: 平顶光闭环反馈自适应优化版（标定映射 ROI）
% 说明：
% 1. 本脚本用于实验电脑，仍会连接相机和 SLM。
% 2. 闭环 ROI 不再依赖形态学找框，也不再只用理想物理公式估计。
% 3. 每次实验前请运行 Pixel_Matching_Probe.m 重新生成 Calibration_Mapping_Data.mat。
% 4. 本脚本会先检查标定质量，再用四点仿射映射划定平顶光 ROI。

clear; clc; close all;

%% 0. 路径与硬件库初始化
solution_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(solution_dir);
addpath(solution_dir);

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

calib_file = fullfile(project_root, './photo/Calibration_Mapping_Data.mat');

if ~libisloaded('SecondDll')
    loadlibrary('SecondDll.dll', 'SecondDll.h');
end

temp_slm_path = fullfile(solution_dir, 'temp_phase_verify.bmp');

%% 1. 强制重置底层硬件
fprintf('--> 正在清理后台进程，准备独占相机...\n');
imaqreset;
pause(1.5);

%% 2. 物理参数与算法网格
N = 1080;
p_slm = 0.0064;         % SLM 像元尺寸，单位 mm
lambda = 0.000532;      % 波长，单位 mm
f = 300;                % 傅里叶透镜焦距，单位 mm
cam_pixel_pitch = 3.45; % 相机像元尺寸，单位 um

slm_full_w = 1920;
slm_full_h = 1080;

D = N * p_slm;
omega = D / 2;
x = linspace(-D/2, D/2, N);
y = linspace(-D/2, D/2, N);
[X, Y] = meshgrid(x, y);

% 焦平面上 1 个算法像素对应的物理间隔，单位 mm。
dp_focal = (lambda * f) / (N * p_slm);

%% 3. 平顶光目标尺寸与算法偏移
target_w_px = 330; % 目标平顶宽度，单位相机像素
target_h_px = 120; % 目标平顶高度，单位相机像素
shift_x = -100;     % 算法平面 X 偏移；正值对应相机 +X
shift_y = -100;     % 算法平面 Y 偏移；正值对应相机 -Y

target_width_um = target_w_px * cam_pixel_pitch;
target_height_um = target_h_px * cam_pixel_pitch;
target_w_algo = round((target_width_um / 1000) / dp_focal);
target_h_algo = round((target_height_um / 1000) / dp_focal);

%% 4. 相机长连接：抓取零级光基准
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
        EXP_MAX = 1000000;
        EXP_MIN = 10;
    end

    src.ExposureTime = 2000;
    start(vid);
    pause(0.5);

    calllib('SecondDll', 'saShowImageFromFilePath', ...
        black_bmp_path, 0, slm_full_w, 0, slm_full_w, slm_full_h, 1);
    pause(0.5);

    for attempt = 1:15
        img_test = getsnapshot(vid);
        current_max = max(max(medfilt2(double(img_test), [3 3])));

        if current_max >= 250
            src.ExposureTime = max(round(src.ExposureTime * 0.6), EXP_MIN);
            pause(0.1);
        elseif current_max < 200
            exposure_gain = min(240 / max(current_max, 1), 2.0);
            src.ExposureTime = min(round(src.ExposureTime * exposure_gain), EXP_MAX);
            pause(0.1);
        else
            break;
        end
    end

    zero_order_img = img_test;
    [~, max_idx] = max(zero_order_img(:));
    [cy, cx] = ind2sub(size(zero_order_img), max_idx);

    fprintf('--> 零级光坐标锁定: (X: %d, Y: %d) | 曝光: %.0f us\n', ...
        cx, cy, src.ExposureTime);
catch ME
    error('相机独占失败: %s', ME.message);
end

%% 5. 计算初始 GS 相位（纯算法迭代 50 次）
fprintf('--> 正在通过 GS 算法计算初始相位图（50 次）...\n');

I_G = exp(-2 * (X.^2 + Y.^2) / omega^2);
Amp_G = sqrt(I_G);

I_T = zeros(N, N);
row_start = round(N/2 - target_h_algo/2) + shift_y + 1;
row_end = round(N/2 + target_h_algo/2) + shift_y;
col_start = round(N/2 - target_w_algo/2) + shift_x + 1;
col_end = round(N/2 + target_w_algo/2) + shift_x;
I_T(row_start:row_end, col_start:col_end) = 1;
I_T = I_T * (sum(I_G(:)) / sum(I_T(:)));
Amp_T = sqrt(I_T);

phi_slm = 2 * pi * rand(N, N);
for k = 1:50
    U_focal = fftshift(fft2(ifftshift(Amp_G .* exp(1i * phi_slm))));
    phi_slm = angle(fftshift(ifft2(ifftshift(Amp_T .* exp(1i * angle(U_focal))))));
end

% 将初始相位写入 SLM 的 1080x1920 画布中心区域。
phi_2pi = mod(phi_slm, 2*pi);
phi_gray_export = zeros(slm_full_h, slm_full_w, 'uint8');
start_x = round((slm_full_w - N) / 2) + 1;
phi_gray_export(:, start_x:start_x+N-1) = uint8((phi_2pi / (2*pi)) * 255);

imwrite(phi_gray_export, temp_slm_path, 'bmp');
calllib('SecondDll', 'saShowImageFromFilePath', ...
    temp_slm_path, 0, slm_full_w, 0, slm_full_w, slm_full_h, 1);
pause(1.0);
phi_slm_initial = phi_slm;
phi_gray_export_initial = phi_gray_export;

%% 6. 相机唤醒：提高曝光准备闭环
fprintf('--> 唤醒相机，自适应提高曝光...\n');

[cam_H, cam_W] = size(zero_order_img);
[CamX, CamY] = meshgrid(1:cam_W, 1:cam_H);
mask_zero = ((CamX - cx).^2 + (CamY - cy).^2) > 150^2;

src.ExposureTime = min(round(src.ExposureTime * 30), EXP_MAX);
pause(0.2);

for attempt = 1:25
    img_test = double(getsnapshot(vid)) .* mask_zero;
    current_max = max(img_test(:));

    if current_max >= 250
        src.ExposureTime = max(round(src.ExposureTime * 0.7), EXP_MIN);
        pause(0.1);
    elseif current_max < 180
        exposure_gain = min(220 / max(current_max, 1), 3.0);
        src.ExposureTime = min(round(src.ExposureTime * exposure_gain), EXP_MAX);
        pause(0.1);
    else
        break;
    end
end

final_exposure = src.ExposureTime;
initial_flattop = getsnapshot(vid);
fprintf('--> 曝光准备完毕: %.0f us。正式进入硬件闭环反馈。\n\n', final_exposure);

%% 7. 硬件在环迭代（Camera-in-the-Loop）
fprintf('================ 开始实验实拍反馈迭代 ================\n');

closed_loop_iters = 50;  % 实验闭环次数
alpha = 0.1;             % 反馈步长，调小可降低震荡风险
feedback_sign = 1;       % 初始保持正反馈符号，若实验发散再单独改为 -1 验证。
A_weight = Amp_T;        % 自适应权重矩阵
Target_Amp_ROI = Amp_T(row_start:row_end, col_start:col_end);
max_sample_oob_ratio = 0.01;
min_sample_mean_intensity = 1.0;
roi_area_sample_count = 7;  % 每个算法 ROI 像素用 7x7 子采样做相机物理面积平均。
diagnostic_save_enabled = closed_loop_iters > 0; % Keep diagnostics on for hardware runs.
diagnostic_max_raw_frames = closed_loop_iters;

% 标定映射 ROI：
% 四点像素匹配结果用于描述当前光路下"算法平面 -> 相机平面"的仿射映射。
% 标定质量不满足阈值时，程序会停止，避免使用错误 ROI 继续闭环。
algo_roi_h = row_end - row_start + 1;
algo_roi_w = col_end - col_start + 1;

[calib_roi, calib_info] = computeCalibrationFlatTopROI( ...
    calib_file, ...
    [cx, cy], ...
    [cam_H, cam_W], ...
    [shift_x, shift_y], ...
    [algo_roi_w, algo_roi_h]);

cam_cx = calib_roi.center_x;
cam_cy = calib_roi.center_y;
cam_r_start = calib_roi.r_start;
cam_r_end = calib_roi.r_end;
cam_c_start = calib_roi.c_start;
cam_c_end = calib_roi.c_end;
cam_roi_h = calib_roi.height;
cam_roi_w = calib_roi.width;
target_mask_cam = calib_roi.mask;

quality = calib_info.quality;
fprintf('--> 已根据当前像素匹配标定锁定平顶光 ROI。\n');
fprintf('    标定文件: %s\n', calib_file);
fprintf('    最大/均方重投影误差: %.3f / %.3f px\n', ...
    quality.max_reprojection_error_px, quality.rms_reprojection_error_px);
fprintf('    尺度误差: %.3f %% | 旋转角: %.3f deg\n', ...
    quality.scale_error_pct, quality.rotation_angle_deg);
fprintf('    预测中心: (X: %d, Y: %d)\n', cam_cx, cam_cy);
fprintf('    预测 ROI: rows [%d, %d], cols [%d, %d], size = %d x %d\n', ...
    cam_r_start, cam_r_end, cam_c_start, cam_c_end, cam_roi_w, cam_roi_h);

if calib_roi.is_clipped
    error('标定映射 ROI 被相机边界截断，请重新检查标定、零级光中心或 shift_x/shift_y 设置。');
end

diag = struct();
if diagnostic_save_enabled
    diag.run_timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    diag.output_dir = fullfile(solution_dir, 'hardware_diagnostics', ['run_' diag.run_timestamp]);
    if ~exist(diag.output_dir, 'dir')
        mkdir(diag.output_dir);
    end

    raw_frame_count = min(closed_loop_iters, diagnostic_max_raw_frames);
    diag.raw_frame_stride = max(1, ceil(closed_loop_iters / raw_frame_count));
    diag.raw_frame_loop = zeros(1, raw_frame_count);
    diag.raw_frame_stack = zeros(cam_H, cam_W, raw_frame_count, 'uint8');
    diag.raw_frame_write_index = 0;

    diag.sampled_roi_raw_stack = zeros(algo_roi_h, algo_roi_w, closed_loop_iters);
    diag.sampled_roi_smooth_stack = zeros(algo_roi_h, algo_roi_w, closed_loop_iters);
    diag.point_sampled_roi_raw_stack = zeros(algo_roi_h, algo_roi_w, closed_loop_iters);
    diag.point_sampled_roi_smooth_stack = zeros(algo_roi_h, algo_roi_w, closed_loop_iters);
    diag.point_normalized_amp_roi_stack = zeros(algo_roi_h, algo_roi_w, closed_loop_iters);
    diag.point_err_rms = nan(1, closed_loop_iters);
    diag.normalized_amp_roi_stack = zeros(algo_roi_h, algo_roi_w, closed_loop_iters);
    diag.weight_roi_stack = zeros(algo_roi_h, algo_roi_w, closed_loop_iters);
    diag.err_rms = nan(1, closed_loop_iters);
    diag.sample_mean = nan(1, closed_loop_iters);
    diag.sample_oob_ratio = nan(1, closed_loop_iters);
    diag.sample_valid_ratio = nan(1, closed_loop_iters);
    diag.sample_x_range = nan(closed_loop_iters, 2);
    diag.sample_y_range = nan(closed_loop_iters, 2);
    diag.exposure_us = nan(1, closed_loop_iters);
end

for loop = 1:closed_loop_iters
    % 1. 获取当前实拍光场。
    img_raw = double(getsnapshot(vid));
    if diagnostic_save_enabled && ...
            (loop == 1 || loop == closed_loop_iters || mod(loop - 1, diag.raw_frame_stride) == 0)
        if diag.raw_frame_write_index < size(diag.raw_frame_stack, 3)
            diag.raw_frame_write_index = diag.raw_frame_write_index + 1;
            diag.raw_frame_loop(diag.raw_frame_write_index) = loop;
            diag.raw_frame_stack(:, :, diag.raw_frame_write_index) = ...
                uint8(max(min(round(img_raw), 255), 0));
        end
    end

    % 2. 通过标定仿射面积平均采样，让误差矩阵对应相机上的物理小区域。
    [cam_sampled_raw, sample_info] = sampleCalibrationROIFromCamera( ...
        img_raw, calib_file, [cx, cy], ...
        row_start, row_end, col_start, col_end, N, roi_area_sample_count);
    cam_resized = imgaussfilt(cam_sampled_raw, 0.6);

    if diagnostic_save_enabled
        [cam_point_sampled_raw, point_sample_info] = sampleCalibrationROIFromCamera( ...
            img_raw, calib_file, [cx, cy], ...
            row_start, row_end, col_start, col_end, N, 1);
        cam_point_resized = imgaussfilt(cam_point_sampled_raw, 0.6);
    end

    if sample_info.out_of_bounds_ratio > max_sample_oob_ratio
        error(['标定反向采样越界比例过高：%.3f%% > %.3f%%。' ...
               '请重新标定，或检查零级光中心/shift_x/shift_y 偏移。'], ...
            100 * sample_info.out_of_bounds_ratio, 100 * max_sample_oob_ratio);
    end

    sample_mean = mean(max(cam_resized(:), 0));
    if sample_mean < min_sample_mean_intensity
        error(['标定反向采样能量过低：mean %.3f < %.3f。' ...
               '请重新标定，或检查光路对准/偏移设置。'], ...
            sample_mean, min_sample_mean_intensity);
    end

    % 3. 提取实拍振幅并按目标 ROI 总振幅归一化。
    A_cam_roi = sqrt(max(cam_resized, 0));
    A_cam_sum = sum(A_cam_roi(:));
    if A_cam_sum <= eps
        error('标定反向采样得到零振幅，请重新标定或检查偏移设置。');
    end
    A_cam_norm = A_cam_roi * (sum(Target_Amp_ROI(:)) / A_cam_sum);

    % 6. 计算 RMS 误差，用于监控闭环收敛情况。
    err_rms = sqrt(mean((Target_Amp_ROI(:) - A_cam_norm(:)).^2));
    if diagnostic_save_enabled
        A_point_roi = sqrt(max(cam_point_resized, 0));
        A_point_sum = sum(A_point_roi(:));
        if A_point_sum > eps
            A_point_norm = A_point_roi * (sum(Target_Amp_ROI(:)) / A_point_sum);
            point_err_rms = sqrt(mean((Target_Amp_ROI(:) - A_point_norm(:)).^2));
        else
            A_point_norm = zeros(size(A_point_roi));
            point_err_rms = NaN;
        end
        fprintf(['  Loop %02d/%02d | area RMS: %.4f | point RMS: %.4f | ' ...
                 'sample oob: %.4f | mean: %.3f\n'], ...
            loop, closed_loop_iters, err_rms, point_err_rms, ...
            sample_info.out_of_bounds_ratio, sample_mean);
    else
        fprintf('  Loop %02d/%02d | area RMS: %.4f | sample oob: %.4f | mean: %.3f\n', ...
            loop, closed_loop_iters, err_rms, sample_info.out_of_bounds_ratio, sample_mean);
    end

    if diagnostic_save_enabled
        diag.sampled_roi_raw_stack(:, :, loop) = cam_sampled_raw;
        diag.sampled_roi_smooth_stack(:, :, loop) = cam_resized;
        diag.point_sampled_roi_raw_stack(:, :, loop) = cam_point_sampled_raw;
        diag.point_sampled_roi_smooth_stack(:, :, loop) = cam_point_resized;
        diag.point_normalized_amp_roi_stack(:, :, loop) = A_point_norm;
        diag.point_err_rms(loop) = point_err_rms;
        diag.normalized_amp_roi_stack(:, :, loop) = A_cam_norm;
        diag.err_rms(loop) = err_rms;
        diag.sample_mean(loop) = sample_mean;
        diag.sample_oob_ratio(loop) = sample_info.out_of_bounds_ratio;
        diag.sample_valid_ratio(loop) = sample_info.valid_pixel_ratio;
        diag.sample_x_range(loop, :) = sample_info.sample_x_range;
        diag.sample_y_range(loop, :) = sample_info.sample_y_range;
        diag.exposure_us(loop) = src.ExposureTime;
    end

    % 7. 更新 Weighted GS 中的局部目标权重。
    Weight_ROI = A_weight(row_start:row_end, col_start:col_end);
    Weight_ROI = Weight_ROI + feedback_sign * alpha * (Target_Amp_ROI - A_cam_norm);
    Weight_ROI(Weight_ROI < 0) = 0;
    A_weight(row_start:row_end, col_start:col_end) = Weight_ROI;
    if diagnostic_save_enabled
        diag.weight_roi_stack(:, :, loop) = Weight_ROI;
    end

    % 8. 用新权重执行少量内部 GS 迭代。
    for k = 1:5
        U_focal = fftshift(fft2(ifftshift(Amp_G .* exp(1i * phi_slm))));
        phi_slm = angle(fftshift(ifft2(ifftshift(A_weight .* exp(1i * angle(U_focal))))));
    end

    % 9. 下发新相位到 SLM。
    phi_2pi = mod(phi_slm, 2*pi);
    phi_gray_export(:, start_x:start_x+N-1) = uint8((phi_2pi / (2*pi)) * 255);
    imwrite(phi_gray_export, temp_slm_path, 'bmp');
    calllib('SecondDll', 'saShowImageFromFilePath', ...
        temp_slm_path, 0, slm_full_w, 0, slm_full_w, slm_full_h, 1);

    pause(0.6);
end

final_flattop = getsnapshot(vid);
fprintf('================ 闭环优化结束 ================\n');

if diagnostic_save_enabled
    diag.raw_frame_loop = diag.raw_frame_loop(1:diag.raw_frame_write_index);
    diag.raw_frame_stack = diag.raw_frame_stack(:, :, 1:diag.raw_frame_write_index);
    [diag.best_err_rms, diag.best_loop] = min(diag.err_rms);

    diag.zero_order_img = uint8(max(min(round(zero_order_img), 255), 0));
    diag.initial_flattop = uint8(max(min(round(initial_flattop), 255), 0));
    diag.final_flattop = uint8(max(min(round(final_flattop), 255), 0));
    diag.Target_Amp_ROI = Target_Amp_ROI;
    diag.A_weight_final = A_weight;
    diag.A_weight_roi_final = A_weight(row_start:row_end, col_start:col_end);
    diag.phi_slm_initial = phi_slm_initial;
    diag.phi_slm_final = phi_slm;
    diag.phi_gray_export_initial = phi_gray_export_initial;
    diag.phi_gray_export_final = phi_gray_export;
    diag.calib_roi = calib_roi;
    diag.calib_info = calib_info;
    calib_snapshot = load(calib_file);
    diag.calib_data = calib_snapshot.calib_data;
    diag.params = struct( ...
        'N', N, ...
        'p_slm', p_slm, ...
        'lambda', lambda, ...
        'f', f, ...
        'cam_pixel_pitch', cam_pixel_pitch, ...
        'target_w_px', target_w_px, ...
        'target_h_px', target_h_px, ...
        'target_w_algo', target_w_algo, ...
        'target_h_algo', target_h_algo, ...
        'shift_x', shift_x, ...
        'shift_y', shift_y, ...
        'row_start', row_start, ...
        'row_end', row_end, ...
        'col_start', col_start, ...
        'col_end', col_end, ...
        'algo_roi_h', algo_roi_h, ...
        'algo_roi_w', algo_roi_w, ...
        'cam_H', cam_H, ...
        'cam_W', cam_W, ...
        'zero_center_xy', [cx, cy], ...
        'final_exposure_us', final_exposure, ...
        'closed_loop_iters', closed_loop_iters, ...
        'alpha', alpha, ...
        'feedback_sign', feedback_sign, ...
        'roi_area_sample_count', roi_area_sample_count, ...
        'max_sample_oob_ratio', max_sample_oob_ratio, ...
        'min_sample_mean_intensity', min_sample_mean_intensity, ...
        'calib_file', calib_file);

    imwrite(diag.zero_order_img, fullfile(diag.output_dir, 'zero_order_img.png'));
    imwrite(diag.initial_flattop, fullfile(diag.output_dir, 'initial_flattop.png'));
    imwrite(diag.final_flattop, fullfile(diag.output_dir, 'final_flattop.png'));
    imwrite(phi_gray_export_initial, fullfile(diag.output_dir, 'phase_initial.bmp'));
    imwrite(phi_gray_export, fullfile(diag.output_dir, 'phase_final.bmp'));
    copyfile(calib_file, fullfile(diag.output_dir, 'Calibration_Mapping_Data.mat'));

    best_raw_index = find(diag.raw_frame_loop == diag.best_loop, 1);
    if ~isempty(best_raw_index)
        imwrite(diag.raw_frame_stack(:, :, best_raw_index), ...
            fullfile(diag.output_dir, sprintf('best_loop_%02d_raw.png', diag.best_loop)));
    end
    imwrite(uint8(255 * mat2gray(diag.sampled_roi_raw_stack(:, :, diag.best_loop))), ...
        fullfile(diag.output_dir, sprintf('best_loop_%02d_area_roi_raw.png', diag.best_loop)));
    imwrite(uint8(255 * mat2gray(diag.sampled_roi_smooth_stack(:, :, diag.best_loop))), ...
        fullfile(diag.output_dir, sprintf('best_loop_%02d_area_roi_smooth.png', diag.best_loop)));
    imwrite(uint8(255 * mat2gray(diag.point_sampled_roi_raw_stack(:, :, diag.best_loop))), ...
        fullfile(diag.output_dir, sprintf('best_loop_%02d_point_roi_raw.png', diag.best_loop)));
    imwrite(uint8(255 * mat2gray(diag.point_sampled_roi_smooth_stack(:, :, diag.best_loop))), ...
        fullfile(diag.output_dir, sprintf('best_loop_%02d_point_roi_smooth.png', diag.best_loop)));

    diag_fig = figure('Name', 'Hardware Loop Diagnostics', ...
        'Position', [80, 80, 1500, 900], ...
        'Color', 'w', ...
        'Visible', 'off');
    tiledlayout(diag_fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(1:closed_loop_iters, diag.err_rms, 'o-', 'LineWidth', 1.4);
    hold on;
    plot(1:closed_loop_iters, diag.point_err_rms, 's--', 'LineWidth', 1.1);
    grid on;
    xlabel('Loop');
    ylabel('RMS error');
    legend({'area average', 'point sample'}, 'Location', 'northeast');
    title(sprintf('Area RMS, best loop %02d = %.4f', diag.best_loop, diag.best_err_rms));

    nexttile;
    plot(1:closed_loop_iters, diag.sample_mean, 'o-', 'LineWidth', 1.4);
    grid on;
    xlabel('Loop');
    ylabel('Mean intensity');
    title('Sampled ROI mean');

    nexttile;
    imshow(diag.initial_flattop, [0, 255]);
    colormap(gca, 'hot');
    hold on;
    rectangle('Position', [cam_c_start, cam_r_start, cam_roi_w, cam_roi_h], ...
        'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '-.');
    title('Initial camera frame');

    nexttile;
    imshow(diag.final_flattop, [0, 255]);
    colormap(gca, 'hot');
    hold on;
    rectangle('Position', [cam_c_start, cam_r_start, cam_roi_w, cam_roi_h], ...
        'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '-.');
    title('Final camera frame');

    nexttile;
    imagesc(diag.sampled_roi_smooth_stack(:, :, diag.best_loop));
    axis image;
    colorbar;
    title('Best area-averaged ROI');

    nexttile;
    imagesc(diag.normalized_amp_roi_stack(:, :, diag.best_loop) - Target_Amp_ROI);
    axis image;
    colorbar;
    title('Best normalized amplitude error');

    exportgraphics(diag_fig, fullfile(diag.output_dir, 'diagnostic_summary.png'), 'Resolution', 150);
    close(diag_fig);

    save(fullfile(diag.output_dir, 'hardware_diagnostics.mat'), 'diag', '-v7.3');
    fprintf('--> 诊断数据已保存: %s\n', diag.output_dir);
end

%% 8. 安全停机与资源释放
if exist('vid', 'var') && isvalid(vid)
    stop(vid);
    delete(vid);
    clear src vid;
end
imaqreset;
fprintf('--> MATLAB 相机句柄已安全释放。\n');

%% 9. 闭环前后对比展示
figure('Name', '纯算法 VS 硬件闭环优化对比', ...
    'Position', [100, 100, 1600, 500], ...
    'Color', 'w');

% 图 1：算法理想预测。
U_sim = fftshift(fft2(ifftshift(Amp_G .* exp(1i * phi_slm))));
subplot(1, 3, 1);
imshow(abs(U_sim).^2, []);
colormap(gca, 'parula');
xlim([N/2 - 200, N/2 + 200]);
ylim([N/2 - 200, N/2 + 200]);
title('1. 算法层面：理论光斑');

% 图 2：闭环前实拍，并画出标定映射 ROI。
subplot(1, 3, 2);
imshow(initial_flattop, [0, 255]);
colormap(gca, 'hot');
hold on;
viscircles([cx, cy], 150, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 0.5);
rectangle('Position', [cam_c_start, cam_r_start, cam_roi_w, cam_roi_h], ...
    'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '-.');
title('2. 实拍：闭环前（标定映射 ROI）');

% 图 3：闭环后实拍，并使用同一个标定映射 ROI。
subplot(1, 3, 3);
imshow(final_flattop, [0, 255]);
colormap(gca, 'hot');
hold on;
viscircles([cx, cy], 150, 'Color', 'w', 'LineStyle', '--', 'LineWidth', 0.5);
rectangle('Position', [cam_c_start, cam_r_start, cam_roi_w, cam_roi_h], ...
    'EdgeColor', 'g', 'LineWidth', 1.5, 'LineStyle', '-.');
title('3. 实拍：闭环反馈 15 次后');
