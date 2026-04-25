% simulate_calibration_flat_top_roi
% 非硬件仿真：验证基于 Calibration_Mapping_Data.mat 的标定 ROI 是否能圈住平顶光。

clear; clc; close all;

solution_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(solution_dir);
addpath(solution_dir);

%% 与 cam_in_loop.m 保持一致的物理参数和目标参数
N = 1080;
p_slm = 0.0064;         % mm
lambda = 0.000532;      % mm
f = 300;                % mm
cam_pixel_pitch = 3.45; % um

target_w_px = 200;
target_h_px = 150;
shift_x = 100;
shift_y = 100;

dp_focal = (lambda * f) / (N * p_slm);
target_width_um = target_w_px * cam_pixel_pitch;
target_height_um = target_h_px * cam_pixel_pitch;
target_w_algo = round((target_width_um / 1000) / dp_focal);
target_h_algo = round((target_height_um / 1000) / dp_focal);

row_start = round(N/2 - target_h_algo/2) + shift_y + 1;
row_end = round(N/2 + target_h_algo/2) + shift_y;
col_start = round(N/2 - target_w_algo/2) + shift_x + 1;
col_end = round(N/2 + target_w_algo/2) + shift_x;
algo_roi_h = row_end - row_start + 1;
algo_roi_w = col_end - col_start + 1;

cam_H = 1080;
cam_W = 1920;

% 仿真用零级中心，只为保证目标落在画幅内；实验时使用实拍 cx/cy。
zero_center_xy = [900, 850];
calib_file = fullfile(project_root, 'Calibration_Mapping_Data.mat');

%% 标定质量检查与 ROI 计算
[calib_roi, calib_info] = computeCalibrationFlatTopROI( ...
    calib_file, ...
    zero_center_xy, ...
    [cam_H, cam_W], ...
    [shift_x, shift_y], ...
    [algo_roi_w, algo_roi_h]);

quality = calib_info.quality;
calib_box = [calib_roi.c_start, calib_roi.r_start, calib_roi.width, calib_roi.height];

fprintf('Calibration file: %s\n', calib_file);
fprintf('Reprojection max/rms = %.6f / %.6f px\n', ...
    quality.max_reprojection_error_px, quality.rms_reprojection_error_px);
fprintf('Scale error = %.6f %% | rotation = %.6f deg\n', ...
    quality.scale_error_pct, quality.rotation_angle_deg);
fprintf('Calibration ROI rows [%d, %d], cols [%d, %d], clipped=%d\n', ...
    calib_roi.r_start, calib_roi.r_end, calib_roi.c_start, calib_roi.c_end, calib_roi.is_clipped);

%% 独立使用标定仿射矩阵生成“真实”模拟平顶光位置
S = load(calib_file);
calib_data = S.calib_data;
tform = affine2d(calib_data.affine_matrix);

half_w = algo_roi_w / 2;
half_h = algo_roi_h / 2;
algo_corners = [
    shift_x - half_w, shift_y - half_h
    shift_x + half_w, shift_y - half_h
    shift_x + half_w, shift_y + half_h
    shift_x - half_w, shift_y + half_h
];

[corner_dx, corner_dy] = transformPointsForward(tform, algo_corners(:, 1), algo_corners(:, 2));
truth_corners = [zero_center_xy(1) + corner_dx, zero_center_xy(2) + corner_dy];

truth_c_start = max(1, floor(min(truth_corners(:, 1))));
truth_c_end = min(cam_W, ceil(max(truth_corners(:, 1))));
truth_r_start = max(1, floor(min(truth_corners(:, 2))));
truth_r_end = min(cam_H, ceil(max(truth_corners(:, 2))));
truth_box = [
    truth_c_start, ...
    truth_r_start, ...
    truth_c_end - truth_c_start + 1, ...
    truth_r_end - truth_r_start + 1];

calib_iou = localBoxIoU(calib_box, truth_box);
fprintf('Calibration ROI / synthetic truth IoU = %.6f\n', calib_iou);

%% 物理公式 ROI 仅作为对照，不作为闭环依据
scale_ratio = (lambda * f) / (N * p_slm * (cam_pixel_pitch / 1000));
physics_cx = round(zero_center_xy(1) + shift_x * scale_ratio);
physics_cy = round(zero_center_xy(2) - shift_y * scale_ratio);
physics_w = max(1, round(algo_roi_w * scale_ratio));
physics_h = max(1, round(algo_roi_h * scale_ratio));
physics_c_start = max(1, physics_cx - round(physics_w / 2));
physics_c_end = min(cam_W, physics_c_start + physics_w - 1);
physics_r_start = max(1, physics_cy - round(physics_h / 2));
physics_r_end = min(cam_H, physics_r_start + physics_h - 1);
physics_box = [
    physics_c_start, ...
    physics_r_start, ...
    physics_c_end - physics_c_start + 1, ...
    physics_r_end - physics_r_start + 1];

physics_iou = localBoxIoU(calib_box, physics_box);
fprintf('Physics formula ROI / calibration ROI IoU = %.6f\n', physics_iou);

%% 生成合成图，并用同一反向采样函数验证像素对应关系
[CamX, CamY] = meshgrid(1:cam_W, 1:cam_H);
rng(11);
sim_img = max(6 + 2 * randn(cam_H, cam_W), 0);
zero_spot = 120 * exp(-((CamX-zero_center_xy(1)).^2 + ...
                         (CamY-zero_center_xy(2)).^2) / (2 * 18^2));
sim_img = sim_img + zero_spot;

inv_tform = invert(tform);
cam_dx_grid = CamX - zero_center_xy(1);
cam_dy_grid = CamY - zero_center_xy(2);
[algo_x_from_cam, algo_y_from_cam] = transformPointsForward(inv_tform, cam_dx_grid, cam_dy_grid);

[algo_col_grid, algo_row_grid] = meshgrid(col_start:col_end, row_start:row_end);
algo_center = (N + 1) / 2;
algo_x_grid = algo_col_grid - algo_center;
algo_y_grid = algo_row_grid - algo_center;

% Use an affine intensity field; bilinear camera sampling should recover it
% to numerical precision when the coordinate mapping is correct.
truth_field_cam = 165 + 0.07 * algo_x_from_cam - 0.05 * algo_y_from_cam;
truth_roi_expected = 165 + 0.07 * algo_x_grid - 0.05 * algo_y_grid;

inside_algo_roi = algo_x_from_cam >= min(algo_x_grid(:)) - 2 & ...
                  algo_x_from_cam <= max(algo_x_grid(:)) + 2 & ...
                  algo_y_from_cam >= min(algo_y_grid(:)) - 2 & ...
                  algo_y_from_cam <= max(algo_y_grid(:)) + 2;
sim_img(inside_algo_roi) = truth_field_cam(inside_algo_roi);
sim_img = min(max(sim_img, 0), 255);

[sampled_roi, sample_info] = sampleCalibrationROIFromCamera( ...
    sim_img, calib_file, zero_center_xy, ...
    row_start, row_end, col_start, col_end, N);
sample_error_rms = sqrt(mean((sampled_roi(:) - truth_roi_expected(:)).^2));

fprintf('Reverse sample out-of-bounds ratio = %.8f\n', sample_info.out_of_bounds_ratio);
fprintf('Reverse sample size = %d x %d (expected %d x %d)\n', ...
    size(sampled_roi, 2), size(sampled_roi, 1), algo_roi_w, algo_roi_h);
fprintf('Reverse sample RMS error = %.10f intensity counts\n', sample_error_rms);

fig = figure('Name', 'Calibration ROI Simulation', 'Color', 'w', 'Visible', 'off');
tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
imagesc(sim_img);
axis image;
colormap(gca, hot);
colorbar;
hold on;

h_zero = plot(zero_center_xy(1), zero_center_xy(2), ...
    'w+', 'MarkerSize', 10, 'LineWidth', 1.5);
rectangle('Position', truth_box, 'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '-');
rectangle('Position', calib_box, 'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--');
rectangle('Position', physics_box, 'EdgeColor', 'c', 'LineWidth', 1.5, 'LineStyle', ':');

h_truth_key = plot(nan, nan, 'y-', 'LineWidth', 2);
h_calib_key = plot(nan, nan, 'g--', 'LineWidth', 2);
h_physics_key = plot(nan, nan, 'c:', 'LineWidth', 1.5);
legend([h_zero, h_truth_key, h_calib_key, h_physics_key], ...
    {'zero order', 'synthetic truth', 'calibration ROI', 'physics reference'}, ...
    'Location', 'southoutside');
title(sprintf('Camera plane: IoU %.4f, max reproj %.3f px', ...
    calib_iou, quality.max_reprojection_error_px));
hold off;

nexttile;
imagesc(sampled_roi);
axis image;
colormap(gca, parula);
colorbar;
title(sprintf('Reverse sampled algorithm ROI: RMS %.3g', sample_error_rms));

out_png = fullfile(solution_dir, 'simulation_calibration_roi_result.png');
exportgraphics(fig, out_png, 'Resolution', 150);
close(fig);
fprintf('Saved simulation figure: %s\n', out_png);

function iou = localBoxIoU(a, b)
    ax1 = a(1); ay1 = a(2); ax2 = a(1) + a(3) - 1; ay2 = a(2) + a(4) - 1;
    bx1 = b(1); by1 = b(2); bx2 = b(1) + b(3) - 1; by2 = b(2) + b(4) - 1;
    ix1 = max(ax1, bx1); iy1 = max(ay1, by1);
    ix2 = min(ax2, bx2); iy2 = min(ay2, by2);
    iw = max(0, ix2 - ix1 + 1);
    ih = max(0, iy2 - iy1 + 1);
    inter = iw * ih;
    area_a = max(0, a(3)) * max(0, a(4));
    area_b = max(0, b(3)) * max(0, b(4));
    iou = inter / max(area_a + area_b - inter, eps);
end
