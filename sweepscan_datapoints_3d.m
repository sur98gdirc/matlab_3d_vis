function sweepscan_datapoints_3d(h_axes, data, select_thickness, ...
    filename_save, n_passes_before_save, n_passes_after_save, ...
    n_frame_edge1, n_frame_edge2, n_steps, sleep_per_step, ...
    frame_origin, frame_edge1_vect, frame_edge2_vect, frame_move_vect)
%

if n_frame_edge1 < 2
    error('number of points along edge 1 of frame must be 2 or more')
end
if n_frame_edge2 < 2
    error('number of points along edge 2 of frame must be 2 or more')
end
if n_steps < 2
    error('number of steps must be 2 or more')
end

if det([ frame_edge1_vect ; frame_edge2_vect ; frame_move_vect ]) == 0
    error('none of frame borders and movement vectors may be collinear')
end

INIT = struct();
INIT.data               = data;
INIT.select_thickness   = select_thickness;
INIT.filename_save      = filename_save;
INIT.n_frame_edge1      = n_frame_edge1;
INIT.n_frame_edge2      = n_frame_edge2;
INIT.n_steps            = n_steps;
INIT.sleep_per_step     = sleep_per_step;
INIT.frame_origin       = frame_origin;
INIT.frame_edge1_vect   = frame_edge1_vect;
INIT.frame_edge2_vect   = frame_edge2_vect;
INIT.frame_move_vect    = frame_move_vect;

PRECOMP = precompute(INIT);

if isempty(h_axes)
    figure();
    h_axes = axes();
end

CURR = struct();
CURR.h_axes = h_axes;
hold('on')
CURR = create_plot_objs(INIT, PRECOMP, CURR);

CURR.filename_save = '';
while n_passes_before_save
    iterate_plot(INIT, PRECOMP, CURR);
    n_passes_before_save = n_passes_before_save - 1;
end

CURR.filename_save = INIT.filename_save;
iterate_plot(INIT, PRECOMP, CURR);

CURR.filename_save = '';
while n_passes_after_save
    iterate_plot(INIT, PRECOMP, CURR);
    n_passes_after_save = n_passes_after_save - 1;
end

return

function PRECOMP = precompute(INIT)
    PRECOMP = struct();
    
    PRECOMP.n_data_pnts = size(INIT.data, 1);
    PRECOMP.data_idxs = 1:PRECOMP.n_data_pnts;
    
    [ PRECOMP.step_vect_by_idx  , PRECOMP.step_idxs  ] = ...
        make_fractions_vector(INIT.n_steps,       INIT.frame_move_vect );
    [ PRECOMP.edge1_vect_by_idx , PRECOMP.edge1_idxs ] = ...
        make_fractions_vector(INIT.n_frame_edge1, INIT.frame_edge1_vect);
    [ PRECOMP.edge2_vect_by_idx , PRECOMP.edge2_idxs ] = ...
        make_fractions_vector(INIT.n_frame_edge2, INIT.frame_edge2_vect);
        
    PRECOMP.origin_by_step      = repmat(INIT.frame_origin, INIT.n_steps, 1) + ...
                                    PRECOMP.step_vect_by_idx;
    PRECOMP.frames_base     = make_frames(INIT, PRECOMP);
    
    PRECOMP.data_xyz            = INIT.data(:, 1:3);
    PRECOMP.distances       = compute_distances(INIT, PRECOMP);
    
    PRECOMP.fl_use_data_point   = PRECOMP.distances <= INIT.select_thickness;
    PRECOMP.frames_data     = fill_frames(INIT, PRECOMP);
    
    PRECOMP.frames_base_p   = permute(PRECOMP.frames_base, [3, 4, 1, 2]);
    PRECOMP.frames_data_p   = permute(PRECOMP.frames_data, [3, 4, 1, 2]);
return

function CURR = create_plot_objs(INIT, PRECOMP, CURR)
    CURR.step_idx = 1;
    CURR.h_surf_base = surf(...
        PRECOMP.frames_base_p(:, :, CURR.step_idx, 1), ...
        PRECOMP.frames_base_p(:, :, CURR.step_idx, 2), ...
        PRECOMP.frames_base_p(:, :, CURR.step_idx, 3), ...
        'EdgeAlpha', 0.2, ...
        'EdgeColor', [0 0 0], ...
        'FaceColor', 'none');
    CURR.h_surf_data = surf(...
        PRECOMP.frames_data_p(:, :, CURR.step_idx, 1), ...
        PRECOMP.frames_data_p(:, :, CURR.step_idx, 2), ...
        PRECOMP.frames_data_p(:, :, CURR.step_idx, 3), ...
        PRECOMP.frames_data_p(:, :, CURR.step_idx, 4), ...
        'UserData', INIT, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
return

function CURR = iterate_plot(INIT, PRECOMP, CURR)
    for step_idx = PRECOMP.step_idxs
        CURR.step_idx = step_idx;
        set(CURR.h_surf_base, 'XData', PRECOMP.frames_base_p(:, :, CURR.step_idx, 1));
        set(CURR.h_surf_base, 'YData', PRECOMP.frames_base_p(:, :, CURR.step_idx, 2));
        set(CURR.h_surf_base, 'ZData', PRECOMP.frames_base_p(:, :, CURR.step_idx, 3));
        set(CURR.h_surf_data, 'XData', PRECOMP.frames_data_p(:, :, CURR.step_idx, 1));
        set(CURR.h_surf_data, 'YData', PRECOMP.frames_data_p(:, :, CURR.step_idx, 2));
        set(CURR.h_surf_data, 'ZData', PRECOMP.frames_data_p(:, :, CURR.step_idx, 3));
        set(CURR.h_surf_data, 'CData', PRECOMP.frames_data_p(:, :, CURR.step_idx, 4));
        if ~isempty(CURR.filename_save)
            saveas(get(CURR.h_axes, 'parent'), sprintf(CURR.filename_save, CURR.step_idx));
        end
        pause(INIT.sleep_per_step);
    end
return

function [vect_by_idx, idxs, fract] = make_fractions_vector(n, vector)
    idxs = 1:n;
    fract = (idxs.' - 1) ./ (n - 1);
    vect_by_idx = fract * vector;
return

function frames_base = make_frames(INIT, PRECOMP)
    frames_base = nan(  INIT.n_steps       , 3 , ...
                        INIT.n_frame_edge1 , ...
                        INIT.n_frame_edge2 );
    
    for         step_idx  = PRECOMP.step_idxs
        for     edge1_idx = PRECOMP.edge1_idxs
                    edge_base                                                       = ...
                    PRECOMP.origin_by_step   (step_idx                         , :) + ...
                    PRECOMP.edge1_vect_by_idx(             edge1_idx           , :);
            for edge2_idx = PRECOMP.edge2_idxs
                    frames_base              (step_idx, :, edge1_idx, edge2_idx   ) = ...
                    edge_base                                                       + ...
                    PRECOMP.edge2_vect_by_idx(                        edge2_idx, :);
                    
            end
        end
    end
return

function distances = compute_distances(INIT, PRECOMP)
    distances = nan(INIT.n_steps, PRECOMP.n_data_pnts);
    for step_idx  = PRECOMP.step_idxs
            distances(step_idx, :) = ...
                dist_plane_points(PRECOMP.origin_by_step(step_idx, :), ...
                    INIT.frame_edge1_vect, INIT.frame_edge2_vect, PRECOMP.data_xyz).';
    end
return

function frames_data = fill_frames(INIT, PRECOMP)
    frames_data = nan(  INIT.n_steps       , 4 , ...
                        INIT.n_frame_edge1 , ...
                        INIT.n_frame_edge2 );
    
    for             step_idx  = PRECOMP.step_idxs
                        dist = inf(INIT.n_frame_edge1, INIT.n_frame_edge2);
        for         data_idx  = PRECOMP.data_idxs
                        if ~ PRECOMP.fl_use_data_point(step_idx, data_idx)
                            continue
                        end
            for     edge1_idx = PRECOMP.edge1_idxs
                for edge2_idx = PRECOMP.edge2_idxs
                        d = norm(PRECOMP.frames_base(step_idx, :,         edge1_idx, edge2_idx) - ...
                                 PRECOMP.data_xyz   (data_idx, :                              ) );
                        if  dist                    (                     edge1_idx, edge2_idx) > d
                            frames_data             (step_idx, :,         edge1_idx, edge2_idx) = ...
                                INIT.data           (data_idx, :                              );
                            dist                    (                     edge1_idx, edge2_idx) = d;
                        end
                end
            end
        end
    end
return

