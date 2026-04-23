%% 1. 相机配置与抓图 (保持之前的优秀逻辑)
clear; clc; close all;

% --- 物理参数输入区 ---
cam_pixel_pitch = 3.45; % 相机像元大小 (um)

% 设定你想要生成的长方形平顶光的物理尺寸 (um)
% 假设你要一个 1200um x 900um 的长方形 (4:3比例)
target_width_um  = 1200; 
target_height_um = 900;  

% 换算成长方形在相机画面上的真实像素个数
target_w_px = round(target_width_um / cam_pixel_pitch);
target_h_px = round(target_height_um / cam_pixel_pitch);

try
    vid = videoinput("gentl", 1, "Mono8");
    src = getselectedsource(vid);
    triggerconfig(vid, 'manual');
    src.ExposureTime = 2000; 
    start(vid);
    pause(0.5); 
    
    fprintf('--> 正在自动调节相机曝光...\n');
    max_attempts = 20;       
    target_max_val = 240;    
    min_acceptable = 200;    
    
    for attempt = 1:max_attempts
        img_test = getsnapshot(vid);
        img_test_filtered = medfilt2(double(img_test), [3 3]); 
        current_max = max(img_test_filtered(:)); 
        
        if current_max >= 250
            src.ExposureTime = src.ExposureTime * 0.6;
            pause(0.1); 
        elseif current_max < min_acceptable
            ratio = target_max_val / max(current_max, 1);
            ratio = min(ratio, 2.0); 
            src.ExposureTime = src.ExposureTime * ratio;
            pause(0.1);
        else
            fprintf('--> 曝光调节成功！最终曝光: %.0f us\n', src.ExposureTime);
            snapshot1 = img_test; 
            img_filtered = img_test_filtered; 
            break;
        end
        if attempt == max_attempts
            warning('曝光未收敛，采用最后一次抓图。');
            snapshot1 = img_test; img_filtered = img_test_filtered;
        end
    end
    stop(vid); delete(vid); clear src vid;
    
catch ME
    warning('使用模拟高斯光斑。');
    [X, Y] = meshgrid(1:1920, 1:1080);
    snapshot1 = uint8(240 * exp(-((X-850).^2 + (Y-500).^2)/(2*20^2)));
    img_filtered = double(snapshot1);
end

%% 2. 焦点定位与 FWHM 计算
[max_val, max_idx] = max(img_filtered(:));
[cy, cx] = ind2sub(size(img_filtered), max_idx);

profile_x = img_filtered(cy, :);
profile_y = img_filtered(:, cx);
half_max = max_val / 2;

idx_x = find(profile_x >= half_max); FWHM_x = idx_x(end) - idx_x(1) + 1;
idx_y = find(profile_y >= half_max); FWHM_y = idx_y(end) - idx_y(1) + 1;
FWHM_mean = round((FWHM_x + FWHM_y) / 2);

fprintf('--> 焦点中心: (X: %d, Y: %d)\n', cx, cy);
fprintf('--> 目标长方形平顶光尺寸: %d x %d 像素\n', target_w_px, target_h_px);

%% 3. 自动划定正方形 ROI (必须包容目标长方形)
% 正方形 ROI 的边长不仅要包住焦点(艾里斑)，还必须包住你想要的长方形目标！
% 策略：取 FWHM的10倍 与 长方形最大边长的1.5倍 之间的最大值
roi_size = max([10 * FWHM_mean, target_w_px * 1.5, target_h_px * 1.5]); 
roi_size = floor(roi_size / 2) * 2; % 保证是偶数
half_roi = roi_size / 2;

% 计算正方形 ROI 边界
x_start = max(1, cx - half_roi);
x_end   = min(size(img_filtered, 2), cx + half_roi - 1);
y_start = max(1, cy - half_roi);
y_end   = min(size(img_filtered, 1), cy + half_roi - 1);

% 提取正方形 ROI 矩阵
img_roi = img_filtered(y_start:y_end, x_start:x_end);

%% 4. 结果可视化分析
figure('Name', 'ROI与目标长方形映射', 'Position', [100, 100, 1200, 500], 'Color', 'w');

% -- 子图 1: 全画幅视图 --
subplot(1, 2, 1);
imshow(snapshot1, []); hold on;
plot(cx, cy, 'r+', 'MarkerSize', 8, 'LineWidth', 1.5); % 光心

% 画绿色虚线框 (正方形 ROI，给算法用的画布)
rectangle('Position', [x_start, y_start, roi_size, roi_size], ...
          'EdgeColor', 'g', 'LineWidth', 2, 'LineStyle', '--');

% 画黄色实线框 (你的长方形平顶光目标)
% 计算长方形在全画幅中的起止点
rect_x = cx - target_w_px/2;
rect_y = cy - target_h_px/2;
rectangle('Position', [rect_x, rect_y, target_w_px, target_h_px], ...
          'EdgeColor', 'y', 'LineWidth', 2);
title('原始全画幅 (绿:算法正方形ROI, 黄:最终平顶光目标)');
legend('初始焦点', '正方形计算区域', '长方形平顶光目标', 'Location', 'southoutside');
hold off;

% -- 子图 2: 提取出的正方形 ROI 视图 (算法真正看到的区域) --
subplot(1, 2, 2);
imshow(uint8(img_roi), []); hold on;
% 在 ROI 坐标系下，中心点就在画面正中央
roi_cx = size(img_roi, 2) / 2;
roi_cy = size(img_roi, 1) / 2;

% 在这个局部坐标系下画出长方形目标
roi_rect_x = roi_cx - target_w_px/2;
roi_rect_y = roi_cy - target_h_px/2;
rectangle('Position', [roi_rect_x, roi_rect_y, target_w_px, target_h_px], ...
          'EdgeColor', 'y', 'LineWidth', 2);
title(sprintf('提取的算法画布 (%d x %d 像素)', roi_size, roi_size));
hold off;