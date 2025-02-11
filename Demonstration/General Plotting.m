%% Deformation Plotting
load("Element_Properties.mat");
load("Real_Deformed.mat");

El_ID=[Element_Start,Element_End];

m=10; % deformation magnification
Nd1=Node_Coor;
Nd1(:,3)=Node_Coor(:,3)-m*(Node_Deformed_1(:,3)-Node_Coor(:,3));
Nd2=Node_Coor;
Nd2(:,3)=Node_Coor(:,3)-m*(Node_Deformed_2(:,3)-Node_Coor(:,3));

G = graph(El_ID(:,1), El_ID(:,2));
axis manual;
view(32,30); % Set 3D view
axis equal;
grid on;
% Plot the graph in 3D with correct node coordinates
h = plot(G, 'XData', Node_Deformed_1(:,1), 'YData', Node_Deformed_1(:,2), 'ZData', Nd1(:,3)); h.NodeLabel = {}; 
h.NodeColor = 'r'; % Set nodes to red
h.EdgeColor = 'b'; % Set edges to blue
h.MarkerSize = 12;
h.LineWidth = 3;
axis manual;
view(31,30);
axis equal;
grid on;

xlim([-2, 12]); ylim([-2, 10]); zlim([-5, 8]);
ax = gca;
ax.Box = 'on';
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
fig = gcf; % Get the current figure handle

% Set the figure background color to white
set(fig, 'Color', 'white');

% Find all text and axes in the figure
all_axes = findall(fig, 'Type', 'axes');
all_text = findall(fig, 'Type', 'text');

% Set font properties for axes
for ax = all_axes'
    set(ax, 'FontName', 'Times', 'FontSize', 30);
end

%% Results Plotting
%% 1. Accuracy
figure;
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
Acc=reshape(Par_Disp,[],1);
Accuracy=[prctile(Acc,25)-median(Acc),median(Acc),prctile(Acc,95)-median(Acc),mean(Acc)]*1000;
Acc=reshape(Gen_Disp,[],1);
Accuracy(2,:)=[prctile(Acc,25)-median(Acc),median(Acc),prctile(Acc,95)-median(Acc),mean(Acc)]*1000;
Acc=reshape(SA_Disp,[],1);
Accuracy(3,:)=[prctile(Acc,25)-median(Acc),median(Acc),prctile(Acc,95)-median(Acc),mean(Acc)]*1000;
categories={'PS','GA','SA'};
Cat = categorical(categories, categories, 'Ordinal', true);
colors = lines(numel(categories));  % Use built-in 'lines' colormap
markers = {'o', 's', 'd'};
hold on;  % Keep multiple plots on the same figure

% Loop through each category to plot separately
for i = 1:numel(categories)
    errorbar(Cat(i),Accuracy(i,2)',Accuracy(i,1)',Accuracy(i,3)','LineWidth', 3, ...
        'Marker', markers{i}, 'MarkerSize', 20, ...
        'MarkerFaceColor', colors(i, :), ...   % Filled marker color
        'MarkerEdgeColor', 'k', ...            % Black border
        'Color', colors(i, :));                % Line color same as marker

end

hold off;

xlabel('Metaheuristic Optimization Method');
ylabel('RMSE (mm)');
grid on;

% Ensure correct ordering
set(gca, 'XTickLabel', categories,'FontSize', 20, 'FontName', 'Times New Roman');
set(gcf, 'Color', 'w');

%% 2. Constraint Violation
nexttile;
Ac=(1-Par_Accuracy(:,5))*100;
Comp=[min(Ac)-median(Ac),median(Ac),max(Ac)-median(Ac)];
Ac=(1-Gen_Accuracy(:,5))*100;
Comp(2,:)=[min(Ac)-median(Ac),median(Ac),max(Ac)-median(Ac)];
Ac=(1-SA_Accuracy(:,5))*100;
Comp(3,:)=[min(Ac)-median(Ac),median(Ac),max(Ac)-median(Ac)];
categories={'PS','GA','SA'};
Cat = categorical(categories, categories, 'Ordinal', true);
colors = lines(numel(categories));  % Use built-in 'lines' colormap
markers = {'o', 's', 'd'};
hold on;  % Keep multiple plots on the same figure

% Loop through each category to plot separately
for i = 1:numel(categories)
    errorbar(Cat(i),Comp(i,2)',Comp(i,1)',Comp(i,3)','LineWidth', 3, ...
        'Marker', markers{i}, 'MarkerSize', 20, ...
        'MarkerFaceColor', colors(i, :), ...   % Filled marker color
        'MarkerEdgeColor', 'k', ...            % Black border
        'Color', colors(i, :));                % Line color same as marker
end

hold off;

xlabel('Metaheuristic Optimization Method');
ylabel('Constraint Violation (%)');
grid on;

% Ensure correct ordering
set(gca, 'XTickLabel', categories,'FontSize', 20, 'FontName', 'Times New Roman');
set(gcf, 'Color', 'w');

%% 3. Plotting the Nodes for Each of the AI-based Algorithms
figure;
t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% 3.1. PS
G = graph(Element_Start,Element_End);
nexttile;
% Create a plot for the graph
h = plot(G, 'XData', Node_Coor(:,1), 'YData', Node_Coor(:,2),'ZData', Node_Coor(:,3));
axis equal
view(31,30);
[b,r]=min(Par_Accuracy(:,4));
% Normalize color data (optional)
C1=Par_Disp(r,:)';
C = (C1 - min(C1)) / (max(C1) - min(C1));  % Normalize C to [0, 1] range
S1=[1;Par_Z_Score(r,:)'];
% Normalize size data (optional)
min_size = 6;
max_size = 18;
S = min_size + ((S1 - min(S1)) * (max_size - min_size)) / (max(S1) - min(S1));  % Scale S
% Set node colors and sizes
h.EdgeColor = 'k';
h.NodeCData = C1;  % Apply color to nodes
h.NodeColor = 'flat';  % Ensure that color data is used
h.MarkerSize = S;  % Apply size to nodes
h.LineWidth = 1;
h.NodeLabel = {};
% Optionally, set a colormap (e.g., 'jet')
colormap(jet);
% colorbar;
% Add colorbar for reference
clim([min(C1), max(C1)]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
% title('PS Node Coordinates');
set(gca,'FontSize', 30, 'FontName', 'Times New Roman');
set(gcf, 'Color', 'w');

% 3.2. GA
G = graph(Element_Start,Element_End);
nexttile;
% Create a plot for the graph
h = plot(G, 'XData', Node_Coor(:,1), 'YData', Node_Coor(:,2),'ZData', Node_Coor(:,3));
axis equal
view(31,30);
[b,r]=min(Gen_Accuracy(:,4));
% Normalize color data (optional)
C1=Gen_Disp(r,:)';
C = (C1 - min(C1)) / (max(C1) - min(C1));  % Normalize C to [0, 1] range
S1=[1;Gen_Z_Score(r,:)'];
% Normalize size data (optional)
min_size = 6;
max_size = 18;
S = min_size + ((S1 - min(S1)) * (max_size - min_size)) / (max(S1) - min(S1));  % Scale S
% Set node colors and sizes
h.EdgeColor = 'k';
h.NodeCData = C1;  % Apply color to nodes
h.NodeColor = 'flat';  % Ensure that color data is used
h.MarkerSize = S;  % Apply size to nodes
h.LineWidth = 1;
h.NodeLabel = {};
% Optionally, set a colormap (e.g., 'jet')
colormap(jet);
% colorbar;
% Add colorbar for reference
clim([min(C1), max(C1)]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
% title('PS Node Coordinates');
set(gca,'FontSize', 30, 'FontName', 'Times New Roman');
set(gcf, 'Color', 'w');

% 3.3. SA
G = graph(Element_Start,Element_End);
nexttile;
% Create a plot for the graph
h = plot(G, 'XData', Node_Coor(:,1), 'YData', Node_Coor(:,2),'ZData', Node_Coor(:,3));
axis equal
view(31,30);
[b,r]=min(SA_Accuracy(:,4));
% Normalize color data (optional)
C1=SA_Disp(r,:)';
C = (C1 - min(C1)) / (max(C1) - min(C1));  % Normalize C to [0, 1] range
S1=[1;SA_Z_Score(r,:)'];
% Normalize size data (optional)
min_size = 6;
max_size = 18;
S = min_size + ((S1 - min(S1)) * (max_size - min_size)) / (max(S1) - min(S1));  % Scale S
% Set node colors and sizes
h.EdgeColor = 'k';
h.NodeCData = C1;  % Apply color to nodes
h.NodeColor = 'flat';  % Ensure that color data is used
h.MarkerSize = S;  % Apply size to nodes
h.LineWidth = 1;
h.NodeLabel = {};
% Optionally, set a colormap (e.g., 'jet')
colormap(jet);

% Add colorbar for reference
% colorbar;
clim([min(C1), max(C1)]);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
% title('PS Node Coordinates');
set(gca,'FontSize', 30, 'FontName', 'Times New Roman');
set(gcf, 'Color', 'w');