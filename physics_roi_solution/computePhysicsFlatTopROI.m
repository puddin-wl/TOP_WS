function [roi, info] = computePhysicsFlatTopROI(zero_center_xy, cam_size_hw, target_center_algo_xy, algo_roi_size_wh, lambda_mm, f_mm, N_algo, p_slm_mm, cam_pixel_pitch_um)
% computePhysicsFlatTopROI Compute the flat-top camera ROI from Fourier optics.
% The calibration MAT file is intentionally not used here. It is only a
% validation artifact for the physical scale relation.

    zero_cx = zero_center_xy(1);
    zero_cy = zero_center_xy(2);
    cam_H = cam_size_hw(1);
    cam_W = cam_size_hw(2);
    shift_x = target_center_algo_xy(1);
    shift_y = target_center_algo_xy(2);
    algo_roi_w = algo_roi_size_wh(1);
    algo_roi_h = algo_roi_size_wh(2);

    cam_pixel_pitch_mm = cam_pixel_pitch_um / 1000;
    dp_focal_mm = (lambda_mm * f_mm) / (N_algo * p_slm_mm);
    scale_ratio = dp_focal_mm / cam_pixel_pitch_mm;

    center_x_float = zero_cx + shift_x * scale_ratio;
    center_y_float = zero_cy - shift_y * scale_ratio;
    cam_cx = round(center_x_float);
    cam_cy = round(center_y_float);

    requested_w = max(1, round(algo_roi_w * scale_ratio));
    requested_h = max(1, round(algo_roi_h * scale_ratio));

    c_start_raw = cam_cx - round(requested_w / 2);
    c_end_raw = c_start_raw + requested_w - 1;
    r_start_raw = cam_cy - round(requested_h / 2);
    r_end_raw = r_start_raw + requested_h - 1;

    c_start = max(1, c_start_raw);
    c_end = min(cam_W, c_end_raw);
    r_start = max(1, r_start_raw);
    r_end = min(cam_H, r_end_raw);

    if c_end < c_start || r_end < r_start
        error('Physics ROI is outside the camera frame. center=(%.2f, %.2f), requested=%dx%d, camera=%dx%d.', ...
            center_x_float, center_y_float, requested_w, requested_h, cam_W, cam_H);
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
    roi.requested_width = requested_w;
    roi.requested_height = requested_h;
    roi.raw_bounds = [r_start_raw, r_end_raw, c_start_raw, c_end_raw];
    roi.is_clipped = r_start ~= r_start_raw || r_end ~= r_end_raw || ...
                     c_start ~= c_start_raw || c_end ~= c_end_raw;
    roi.mask = false(cam_H, cam_W);
    roi.mask(r_start:r_end, c_start:c_end) = true;

    info.dp_focal_mm = dp_focal_mm;
    info.cam_pixel_pitch_mm = cam_pixel_pitch_mm;
    info.scale_ratio = scale_ratio;
    info.formula = 'scale_ratio = (lambda * f) / (N * p_slm * p_cam)';
    info.sign_convention = 'camera_x = zero_x + shift_x * scale; camera_y = zero_y - shift_y * scale';
end
