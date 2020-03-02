clear all; close all;
warning('off','all')

shape = 'triangle';
MAX_NEIGHBORS = 20;
CIRCLE_SIDES = 30;
THETA = linspace(0, 2*pi, CIRCLE_SIDES);
PERIM_THRESHOLD = 0.00001;

polygon = fopen(shape, 'r');
formatSpec = '%f %f';
sizeA = [2 Inf];
P = fscanf(polygon, formatSpec, sizeA);
fclose(polygon);
% close off the polygon
plgn_len = length(P);

med_axis_nodes = fopen(strcat(shape, "_med_axis_nodes"), 'r');
formatSpec = '%d %f %f';
sizeA = [3 Inf];
V_id_pt = fscanf(med_axis_nodes, formatSpec, sizeA);
fclose(med_axis_nodes);
V_ids = V_id_pt(1,:);
V = V_id_pt(2:3,:);

% display(V_ids);
% display(V);

med_axis_edges = fopen(strcat(shape,'_med_axis_edges'), 'r');
formatSpec = '%d %d';
sizeA = [2 Inf];
E = fscanf(med_axis_edges, formatSpec, sizeA);
fclose(med_axis_edges);

% create line segments to visualize medial axis
num_edges = length(E);
mdax_xs = zeros(2,num_edges);
mdax_ys = zeros(2,num_edges);
for i = 1:num_edges
    id_pt1 = find(V_ids==E(1,i));
    id_pt2 = find(V_ids==E(2,i));
    mdax_xs(:,i) = [V(1, id_pt1); V(1, id_pt2)];
    mdax_ys(:,i) = [V(2, id_pt1); V(2, id_pt2)];
end

% visualize the polygon line segments
shifted_P = [P(:,2:plgn_len), P(:,1)];
plgn_xs = [P(1,:); shifted_P(1,:)];
plgn_ys = [P(2,:); shifted_P(2,:)];

sampled_nodes = fopen(strcat(shape, "_sampled_axis_nodes"), 'r');
formatSpec = '%d %f %f %f';
sizeA = [4 Inf];
V_id_pt = fscanf(sampled_nodes, formatSpec, sizeA);
fclose(sampled_nodes);
V_ids = uint64(V_id_pt(1,:));
V = V_id_pt(2:3,:);
R = V_id_pt(4,:);
N = length(V_ids);

sampled_edges = fopen(strcat(shape,'_sampled_axis_edges'), 'r');
formatSpec = '%d %d';
sizeA = [2 Inf];
E = fscanf(sampled_edges, formatSpec, sizeA);
fclose(sampled_edges);

% create vertex index map
id_map = containers.Map(V_ids, 1:N);

% create neighbor matrix
neighbors = zeros(N, MAX_NEIGHBORS);

% find peaks of radius function by bfs
% also create edge map
radius_peaks = zeros(1,N);
peak_count = 0;
E_full = [E, [E(2,:);E(1,:)]];
for id = V_ids
    indx = id_map(id);
    radius = R(indx);
    node_neighbors = E_full(2, E_full(1,:)==id);
    
    % insert into neighbor matrix
    num_neighbors = length(node_neighbors);
    neighbors(indx, 1:num_neighbors) = node_neighbors;
    
    is_peak = true;
    for neighbor = node_neighbors
        neighbor_idx = id_map(neighbor);
        neighbor_radius = R(neighbor_idx);
        if neighbor_radius >= radius
            is_peak = false;
            break
        end
    end
    if is_peak
        radius_peaks(indx) = 1;
    end
end

% disp(V_ids(radius_peaks==1));
% disp(R(id_map(1)));
% neighbors = E_full(2, E_full(1,:)==1);
% disp(neighbors);
% for j = neighbors
%     disp(j)
%     disp(R(id_map(j)));
% end

% vertices included in approximation
included_count = 0;
included = zeros(1, N);

% store perimeter/area of Cheeger sets
perims = zeros(1, N);
areas = zeros(1,N);

% candidates for next addition
horizon = zeros(1,N);

% find starting point of iteration
[max_rad, max_indx] = max(R);
included(max_indx) = 1;
included_count = included_count + 1;

for k = neighbors(max_indx,:)
    if k == 0
        break
    end
    horizon(id_map(k)) = 1;
end
for k = find(radius_peaks)
    horizon(k) = 1;
end

% set up all the circles
% take unions at each step
% kind of inefficient, but sometimes MATLAB
% messes up and doesn't compute everything

circles = repelem(polyshape(), N);
for k = V_ids
    k_indx = id_map(k);
    center = V(:,k_indx);
    radius = R(k_indx);
    xys = [cos(THETA); sin(THETA)]*radius + center;
    circles(k_indx) = polyshape(xys(1,:), xys(2,:));
end

% set up first polygon
current_polyshape = union(circles(included>0));
%current_polyshape = union(circles);
current_area = perimeter(current_polyshape);
current_perim = area(current_polyshape);

f = figure('visible', 'off');

while true
    if mod(included_count, 15) == 1
        hold on
        axis equal
        plot(mdax_xs, mdax_ys, 'r')
        plot(plgn_xs, plgn_ys, 'black')
        plot(union(circles(included>0)));
        saveas(f,sprintf('results_%s_%d_iters.png', shape, included_count-1));
        hold off
        clf(f);
    end
    max_ratio = -Inf;
    max_indx = -Inf;
    % max_polyshape = polyshape();
    max_area = -Inf;
    max_perim = -Inf;
    frontier = find(horizon);
    % if no answers have valid ratios, accept any nan ratio
    nan_index_store = -1;
    if isempty(frontier)
        disp("Empty Frontier")
        break
    end
    for p = frontier
        %close_perim_flag = false;
        %iterate through neighbors
        included(p) = 1;
        possible_polyshape = union(circles(included>0));
        possible_perim = perimeter(possible_polyshape);
        possible_area = area(possible_polyshape);
        diff_ratio = (possible_area - current_area)/(possible_perim - current_perim);
        included(p) = 0;
        if isnan(diff_ratio) && nan_index_store < 0
            nan_index_store = p;
            max_area = possible_area;
            max_perim = possible_perim;
        elseif diff_ratio > max_ratio
            max_ratio = diff_ratio;
            max_indx = p;
            max_area = possible_area;
            max_perim = possible_perim;
        end
    end
    
    % do all necessary area and perim updates
    if max_indx >= 0
        included(max_indx) = 1;
        included_count = included_count + 1;
        areas(included_count) = max_area;
        perims(included_count) = max_perim;
    elseif nan_index_store >= 0
        included(max_indx) = 1;
        included_count = included_count + 1;
        areas(included_count) = max_area;
        perims(included_count) = max_perim;
    else
        break
    end
    
    % update horizon
    horizon(max_indx) = 0; % take the node out of neighbor set
    for i = neighbors(max_indx,:)
        if i == 0
            continue %change later to break if works
        end
        k = id_map(i);
        if included(k) == 0
            horizon(k) = 1;
        end
    end
end


hold on
axis equal
plot(mdax_xs, mdax_ys, 'r')
plot(plgn_xs, plgn_ys, 'black')
plot(union(circles(included>0)));
saveas(f,sprintf('results_%s_final_iter.png', shape));
hold off
clf(f);


scatter(perims, areas)
saveas(f,sprintf('results_%s_scatter.png', shape));

