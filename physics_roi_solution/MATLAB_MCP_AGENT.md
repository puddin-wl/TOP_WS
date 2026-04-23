# MATLAB MCP Agent Log

- Task: Build a physics-mapping ROI version of cam_in_loop.m and a non-hardware simulation.
- Correction: Calibration_Mapping_Data.mat is used only to validate the physical scale, not to drive the closed-loop ROI.
- Physical relation: scale_ratio = (lambda * f) / (N * p_slm * p_cam).
- Sign convention: camera_x = zero_x + shift_x * scale; camera_y = zero_y - shift_y * scale.
- Hardware script copied into this folder and modified locally; root cam_in_loop.m was restored to the original version.
- Validation should use MATLAB Code Analyzer and simulate_physics_flat_top_roi.m on non-experimental machines.
