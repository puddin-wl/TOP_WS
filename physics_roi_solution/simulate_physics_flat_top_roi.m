% simulate_physics_flat_top_roi
% Non-hardware simulation for validating the physics-based flat-top ROI.
clear; clc; close all;

solution_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(solution_dir);
addpath(solution_dir);

% Same physical and target parameters used by cam_in_loop.m.
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
row_end   = round(N/2 + target_h_algo/2) + shift_y;
col_start = round(N/2 - target_w_algo/2) + shift_x + 1;
col_end   = round(N/2 + target_w_algo/2) + shift_x;
algo_roi_h = row_end - row_start + 1;
algo_roi_w = col_end - col_start + 1;

cam_H = 1080;
cam_W = 1920;
% Chosen only for simulation so the shifted target remains in-frame.
zero_center_xy = [900, 850];

[pred_roi, pred_info] = computePhysicsFlatTopROI( ...
    zero_center_xy, [cam_H, cam_W], [shift_x, shift_y], ...
    [algo_roi_w, algo_roi_h], lambda, f, N, p_slm, cam_pixel_pitch);

% Build the synthetic image from the same physics relation, independently of the ROI helper.
scale_ratio = (lambda * f) / (N * p_slm * (cam_pixel_pitch / 1000));
truth_cx = round(zero_center_xy(1) + shift_x * scale_ratio);
truth_cy = round(zero_center_xy(2) - shift_y * scale_ratio);
truth_w = max(1, round(algo_roi_w * scale_ratio));
truth_h = max(1, round(algo_roi_h * scale_ratio));
truth_c_start = max(1, truth_cx - round(truth_w / 2));
truth_c_end = min(cam_W, truth_c_start + truth_w - 1);
truth_r_start = max(1, truth_cy - round(truth_h / 2));
truth_r_end = min(cam_H, truth_r_start + truth_h - 1);

[CamX, CamY] = meshgrid(1:cam_W, 1:cam_H);
rng(7);
sim_img = 6 + 2 * randn(cam_H, cam_W);
sim_img = max(sim_img, 0);
zero_spot = 120 * exp(-((CamX-zero_center_xy(1)).^2 + (CamY-zero_center_xy(2)).^2) / (2 * 18^2));
sim_img = sim_img + zero_spot;
flat_patch = 165 + 8 * randn(truth_r_end-truth_r_start+1, truth_c_end-truth_c_start+1);
sim_img(truth_r_start:truth_r_end, truth_c_start:truth_c_end) = flat_patch;
sim_img = min(max(sim_img, 0), 255);

pred_box = [pred_roi.c_start, pred_roi.r_start, pred_roi.width, pred_roi.height];
truth_box = [truth_c_start, truth_r_start, truth_c_end-truth_c_start+1, truth_r_end-truth_r_start+1];
iou_physics = localBoxIoU(pred_box, truth_box);

fprintf('Physics scale_ratio = %.9f camera px / algorithm px\n', pred_info.scale_ratio);
fprintf('Algorithm ROI size = %d x %d algorithm px\n', algo_roi_w, algo_roi_h);
fprintf('Predicted camera ROI = %d x %d camera px\n', pred_roi.width, pred_roi.height);
fprintf('Predicted ROI rows [%d, %d], cols [%d, %d], clipped=%d\n', ...
    pred_roi.r_start, pred_roi.r_end, pred_roi.c_start, pred_roi.c_end, pred_roi.is_clipped);
fprintf('Synthetic truth box rows [%d, %d], cols [%d, %d]\n', ...
    truth_r_start, truth_r_end, truth_c_start, truth_c_end);
fprintf('Physics ROI / synthetic truth IoU = %.6f\n', iou_physics);

calib_file = fullfile(project_root, 'Calibration_Mapping_Data.mat');
has_calib = exist(calib_file, 'file') == 2;
calib_box = [];
calib_iou = NaN;
if has_calib
    S = load(calib_file);
    calib_data = S.calib_data;
    fprintf('Calibration actual_scale = %.9f, theory_scale = %.9f, rotation = %.6f deg\n', ...
        calib_data.actual_scale, calib_data.theory_scale, calib_data.rotation_angle_deg);
    tform = affine2d(calib_data.affine_matrix);
    half_w = algo_roi_w / 2;
    half_h = algo_roi_h / 2;
    algo_corners = [shift_x-half_w, shift_y-half_h; shift_x+half_w, shift_y-half_h; ...
                    shift_x+half_w, shift_y+half_h; shift_x-half_w, shift_y+half_h];
    [corner_dx, corner_dy] = transformPointsForward(tform, algo_corners(:,1), algo_corners(:,2));
    calib_corners = [zero_center_xy(1)+corner_dx, zero_center_xy(2)+corner_dy];
    calib_x1 = floor(min(calib_corners(:,1)));
    calib_x2 = ceil(max(calib_corners(:,1)));
    calib_y1 = floor(min(calib_corners(:,2)));
    calib_y2 = ceil(max(calib_corners(:,2)));
    calib_box = [calib_x1, calib_y1, calib_x2-calib_x1+1, calib_y2-calib_y1+1];
    calib_iou = localBoxIoU(pred_box, calib_box);
    fprintf('Calibration diagnostic box / physics ROI IoU = %.6f\n', calib_iou);
end

fig = figure('Name', 'Physics ROI Simulation', 'Color', 'w', 'Visible', 'off');
imagesc(sim_img); axis image; colormap hot; colorbar; hold on;
h_zero = plot(zero_center_xy(1), zero_center_xy(2), 'w+', 'MarkerSize', 10, 'LineWidth', 1.5);
rectangle('Position', truth_box, 'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '-');
rectangle('Position', pred_box, 'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--');
h_truth_key = plot(nan, nan, 'y-', 'LineWidth', 2);
h_physics_key = plot(nan, nan, 'g--', 'LineWidth', 2);
if has_calib
    rectangle('Position', calib_box, 'EdgeColor', 'c', 'LineWidth', 1.5, 'LineStyle', ':');
    h_calib_key = plot(nan, nan, 'c:', 'LineWidth', 1.5);
    legend([h_zero, h_truth_key, h_physics_key, h_calib_key], ...
        {'zero order', 'synthetic truth', 'physics ROI', 'calibration diagnostic'}, ...
        'Location', 'southoutside');
else
    legend([h_zero, h_truth_key, h_physics_key], ...
        {'zero order', 'synthetic truth', 'physics ROI'}, ...
        'Location', 'southoutside');
end
title(sprintf('Physics ROI validation: IoU %.4f, scale %.6f', iou_physics, pred_info.scale_ratio));
hold off;
out_png = fullfile(solution_dir, 'simulation_physics_roi_result.png');
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
