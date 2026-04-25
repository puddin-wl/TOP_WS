function [roi, info] = computeCalibrationFlatTopROI(calib_file, zero_center_xy, cam_size_hw, target_center_algo_xy, algo_roi_size_wh, thresholds)
% computeCalibrationFlatTopROI 使用四点标定仿射映射生成平顶光相机 ROI。
%
% 输入坐标约定：
%   zero_center_xy        当前实验实拍零级光中心 [cx, cy]
%   target_center_algo_xy 算法平面目标中心偏移 [shift_x, shift_y]
%   algo_roi_size_wh      算法目标区域尺寸 [algo_roi_w, algo_roi_h]

    if nargin < 6
        thresholds = [];
    end

    [calib_data, quality] = validateCalibrationMapping(calib_file, thresholds);
    if ~quality.is_valid
        error(['当前像素匹配标定质量不满足闭环要求，请重新运行 Pixel_Matching_Probe.m。' newline ...
               'max reprojection = %.3f px (limit %.3f), scale error = %.3f %% (limit %.3f), rotation = %.3f deg (limit %.3f).'], ...
            quality.max_reprojection_error_px, quality.thresholds.max_reprojection_error_px, ...
            quality.scale_error_pct, quality.thresholds.max_scale_error_pct, ...
            quality.rotation_angle_deg, quality.thresholds.max_rotation_abs_deg);
    end

    zero_cx = zero_center_xy(1);
    zero_cy = zero_center_xy(2);
    cam_H = cam_size_hw(1);
    cam_W = cam_size_hw(2);
    shift_x = target_center_algo_xy(1);
    shift_y = target_center_algo_xy(2);
    algo_roi_w = algo_roi_size_wh(1);
    algo_roi_h = algo_roi_size_wh(2);

    tform = makeAffineTransform(calib_data.affine_matrix);

    [center_dx, center_dy] = transformPointsForward(tform, shift_x, shift_y);
    center_x_float = zero_cx + center_dx;
    center_y_float = zero_cy + center_dy;
    cam_cx = round(center_x_float);
    cam_cy = round(center_y_float);

    half_w = algo_roi_w / 2;
    half_h = algo_roi_h / 2;
    algo_corners = [
        shift_x - half_w, shift_y - half_h
        shift_x + half_w, shift_y - half_h
        shift_x + half_w, shift_y + half_h
        shift_x - half_w, shift_y + half_h
    ];

    [corner_dx, corner_dy] = transformPointsForward(tform, algo_corners(:, 1), algo_corners(:, 2));
    cam_corners = [zero_cx + corner_dx, zero_cy + corner_dy];

    c_start_raw = floor(min(cam_corners(:, 1)));
    c_end_raw = ceil(max(cam_corners(:, 1)));
    r_start_raw = floor(min(cam_corners(:, 2)));
    r_end_raw = ceil(max(cam_corners(:, 2)));

    c_start = max(1, c_start_raw);
    c_end = min(cam_W, c_end_raw);
    r_start = max(1, r_start_raw);
    r_end = min(cam_H, r_end_raw);

    if c_end < c_start || r_end < r_start
        error('标定映射 ROI 完全落在相机画幅之外，请重新检查标定或 shift_x/shift_y。');
    end

    roi.center_x = cam_cx;
    roi.center_y = cam_cy;
    roi.center_x_float = center_x_float;
    roi.center_y_float = center_y_float;
    roi.r_start = r_start;
    roi.r_end = r_end;
    roi.c_start = c_start;
    roi.c_end = c_end;
    roi.width = c_end - c_start + 1;
    roi.height = r_end - r_start + 1;
    roi.raw_bounds = [r_start_raw, r_end_raw, c_start_raw, c_end_raw];
    roi.is_clipped = r_start ~= r_start_raw || r_end ~= r_end_raw || ...
                     c_start ~= c_start_raw || c_end ~= c_end_raw;
    roi.cam_corners = cam_corners;
    roi.mask = false(cam_H, cam_W);
    roi.mask(r_start:r_end, c_start:c_end) = true;

    info.source = 'Calibration_Mapping_Data.mat / calib_data.affine_matrix';
    info.quality = quality;
    info.actual_scale = calib_data.actual_scale;
    info.theory_scale = calib_data.theory_scale;
    info.rotation_angle_deg = calib_data.rotation_angle_deg;
end

function tform = makeAffineTransform(T)
    if isa(T, 'affine2d')
        tform = T;
        return;
    end

    if ~isequal(size(T), [3, 3])
        error('calib_data.affine_matrix 必须是 affine2d 或 3x3 仿射矩阵。');
    end
    tform = affine2d(T);
end
