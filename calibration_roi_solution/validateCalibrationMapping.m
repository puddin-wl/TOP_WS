function [calib_data, report] = validateCalibrationMapping(calib_file, thresholds)
% validateCalibrationMapping 读取并检查四点像素匹配标定质量。
%
% 标定文件必须由 Pixel_Matching_Probe.m 生成，并包含 calib_data 结构。
% 默认严格阈值：
%   最大重投影误差 <= 3 px
%   尺度误差绝对值 <= 5 %
%   旋转角绝对值 <= 5 deg

    if nargin < 2 || isempty(thresholds)
        thresholds = struct();
    end
    thresholds = fillDefaultThresholds(thresholds);

    if ~exist(calib_file, 'file')
        error('未找到标定文件: %s。请先运行 Pixel_Matching_Probe.m。', calib_file);
    end

    S = load(calib_file);
    if ~isfield(S, 'calib_data')
        error('标定文件缺少 calib_data 结构，请重新运行 Pixel_Matching_Probe.m。');
    end
    calib_data = S.calib_data;

    required_fields = { ...
        'algo_pts', ...
        'cam_pts_centered', ...
        'affine_matrix', ...
        'actual_scale', ...
        'theory_scale', ...
        'rotation_angle_deg'};

    for k = 1:numel(required_fields)
        if ~isfield(calib_data, required_fields{k})
            error('calib_data 缺少字段 "%s"，请重新运行 Pixel_Matching_Probe.m。', required_fields{k});
        end
    end

    algo_pts = calib_data.algo_pts;
    cam_pts_centered = calib_data.cam_pts_centered;
    if size(algo_pts, 2) ~= 2 || size(cam_pts_centered, 2) ~= 2 || ...
            size(algo_pts, 1) ~= size(cam_pts_centered, 1)
        error('标定点尺寸不正确：algo_pts 和 cam_pts_centered 必须是同样行数的 Nx2 数组。');
    end

    tform = makeAffineTransform(calib_data.affine_matrix);
    [pred_x, pred_y] = transformPointsForward(tform, algo_pts(:, 1), algo_pts(:, 2));
    pred_pts = [pred_x, pred_y];

    residuals = pred_pts - cam_pts_centered;
    reproj_errors = sqrt(sum(residuals.^2, 2));
    scale_error_pct = (calib_data.actual_scale - calib_data.theory_scale) / ...
        calib_data.theory_scale * 100;

    report.thresholds = thresholds;
    report.predicted_cam_pts_centered = pred_pts;
    report.residuals = residuals;
    report.reprojection_errors_px = reproj_errors;
    report.max_reprojection_error_px = max(reproj_errors);
    report.mean_reprojection_error_px = mean(reproj_errors);
    report.rms_reprojection_error_px = sqrt(mean(reproj_errors.^2));
    report.scale_error_pct = scale_error_pct;
    report.rotation_angle_deg = calib_data.rotation_angle_deg;

    report.pass_reprojection = report.max_reprojection_error_px <= thresholds.max_reprojection_error_px;
    report.pass_scale = abs(report.scale_error_pct) <= thresholds.max_scale_error_pct;
    report.pass_rotation = abs(report.rotation_angle_deg) <= thresholds.max_rotation_abs_deg;
    report.is_valid = report.pass_reprojection && report.pass_scale && report.pass_rotation;
end

function thresholds = fillDefaultThresholds(thresholds)
    if ~isfield(thresholds, 'max_reprojection_error_px')
        thresholds.max_reprojection_error_px = 3;
    end
    if ~isfield(thresholds, 'max_scale_error_pct')
        thresholds.max_scale_error_pct = 5;
    end
    if ~isfield(thresholds, 'max_rotation_abs_deg')
        thresholds.max_rotation_abs_deg = 5;
    end
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
