clear all; close all;

shape = 'sq_dumbll';

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


hold on
plot(mdax_xs, mdax_ys, 'r')
plot(plgn_xs, plgn_ys, 'black')
hold off



