clc
clearvars

%collection date for verification
images= load('rotation_data.mat')

R=images.ref;
C=images.cur;

figure('Name','rotation data reference image','NumberTitle','off');
imshow(R)
figure('Name','rotation data current image','NumberTitle','off');
imshow(C)

%%%%%%%%%%%Step 2: Exclude Region and define ROI %%%%%%%%%%

%%%%delete/supress this for no exclusion-start, follow line where the stop indication%%%%%%%%%
figure;
imshow(R);
title('Select the region to be excluded in the grid');

%create a mask to exclude the undeformed region 
mask_for_exclude = false(size(R, 1), size(R, 2));

% Loop for interactive selection, replace N with the number of regions to exclude
%total 3 or 4 or 5 regeion will be excluded
%double click after completing and enclosed region

for i = 1:2% the second value might be equal to number of excluded area
    h = drawfreehand('Label', ['Region ' num2str(i)]);  
    wait(h); % hold until all selection is completed
    mask = createMask(h); % Create a mask from the selection
    mask_for_exclude = mask_for_exclude | mask; % Combine the current mask with the previous masks
end

%ROI=region of interest
ROI= ~mask_for_exclude; 

%showing region of interest
figure;
imshow(R); 
hold on; 
h = imshow(~ROI);  % Assuming false in roiMask means exclusion
set(h, 'AlphaData', 0.3);  % Making the mask semi-transparent
title('Excluded and ROI on reference image');
%%%%%%%%%%delete/supress this for no exclusion-end%%%%%%%%


%%%%%%%%%%%%%%%%%%%Step 3: Meshgrid & initial variable for Iteration-Plot%%%%%%%%%%%%%%%%%

%low grid_spacing give more accurate results but increase computational cost
grid_spaceing=20; %two grid points distance
[rows,columns] = size(R);
[xGrid, yGrid] = meshgrid(1:grid_spaceing:columns, 1:grid_spaceing:rows);


figure('Name','Reference Image Grid','NumberTitle','off');
ax = axes; % Create axes in the figure
imshow(R, 'Parent', ax); % Display the image in these axes

figure('Name','Moved subset in deformed image','NumberTitle','off');
ay = axes; % Create axes in the figure
imshow(C, 'Parent', ay); % Display the image in these axes


w=10; % Width-horizontal length-cartesian x
h=10; %height-vertical lenght-cartesian y

R_template_size=[w,h];

x_displacement = zeros(size(xGrid));
y_displacement = zeros(size(yGrid));
final_displacement =zeros(size(xGrid));

%%%%%%%%%%%%%%%%%Step 4: normxcorr2: displacement%%%%%%%%%%%%%%%%%%%%%%%%%%

% Iterate through the grid points
%xGrid and y%Grid have the same dimension (2D)
for i = 1:size(xGrid, 1)
    for j = 1:size(xGrid, 2)
        %if ROI(i, j)  % If the grid point is within the ROI
        %if ~mask_for_exclude(round(yGrid(i)), round(xGrid(i)))  
        %if ROI(round(yGrid(i)), round(xGrid(j)))-ok
        %if ROI((yGrid(i)), (xGrid(j)))
            x=xGrid(i,j);%
            y=yGrid(i,j);%

            % Define the template from the undeformed image
            xMin = max(1, x - R_template_size(1)/2);
            xMax = min(columns, x + R_template_size(1)/2 -1);
            yMin = max(1, y - R_template_size(2)/2);
            yMax = min(rows, y + R_template_size(2)/2-1);

            % Ensure the entire template is within the ROI    
         if all(ROI(yMin:yMax, xMin:xMax), 'all')  % This line would be deleted/supress during no exclusion
            R_template = R(yMin:yMax, xMin:xMax);
 
            % Perform cross-correlation and find the peak(c) in each iteration
            c = normxcorr2(R_template, C);
            [ypeak, xpeak] = find(c==max(c(:))); 

            %subset were defined from middle position so as the yoffset
            yoffSet = ypeak-(size(R_template,1)/2)+1;  
            xoffSet = xpeak-(size(R_template,2)/2)+1;  
             
            % Calculate the displacement
            x_displacement(i, j) = xoffSet - x;
            y_displacement(i, j) = yoffSet - y;

            final_displacement(i,j) =sqrt(x_displacement(i,j).^2 +y_displacement(i,j).^2);
           
            %plotting subset in the terminal
          
            %plot of grid of R
            rectangle(ax, 'Position', [x, y, size(R_template,2), size(R_template,1)], 'EdgeColor', 'g','LineWidth', .5);

            %plot of deformed grid position on c 
            rectangle(ay, 'Position', [xoffSet, yoffSet, size(R_template,2), size(R_template,1)], 'EdgeColor', 'r','LineWidth', .5);
   
            %plot of grid and deformed in current coordinate
            rectangle(ay, 'Position', [x, y, size(R_template,2), size(R_template,1)], 'EdgeColor', 'g','LineWidth', 1);

          
        end%% This line would be deleted/supress during no exclusion
    end
end

%%%%%%%%%%%%%%%%%Step 5: Displacement Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%all displacement plot as well as the next strain plot can hanldle the zero displaement or same dsiplacement in each point plot without showing any error. 
%so if this code searh standard deviation of the data before activate countourf function than it can handle almost all validation process without error
%using countourf function of matlab, one of the requisite is that all input array must be in same dimension
% no variation on array value will not defined in countourf

%xdisplacement plot
fprintf('The x_displacement array are:');
X_Displacement= x_displacement

if std(x_displacement(:) ~= 0) 

     %standard Deviation(std) must not be zero for countourf plot 
     %carefully select ROI so that no undeformed region is chosen.
     %check workspace varible of the corresponding date whether the deviation is zero or not
     % at least one value of x_displacement must not be zero for plotting active countourf
     % i.e no displacement - no graph.No error during valiadation process with same image, translation, rotation of the image
     % Manually define contour levels around the range of interest for better visualization

     figure('Name','xdisplacement contour plot','NumberTitle','off');
     levels = linspace(min(x_displacement(:)), max(x_displacement(:)), 10); % see thw workspace variable for x displacement value range. if it is small use low number 
     contourf(xGrid, yGrid, x_displacement, levels, 'LineColor', 'r');
     colorbar; % Adds a color bar to indicate the displacement values
     title('X Displacement Contour Plot');
     %colormap('parula'); % This can make small differences more visible colormap such as jet, hot, parula.
     xlabel('X-pixel');
     ylabel('Y -pixel');
     axis equal; % original size
     set(gca, 'YDir', 'reverse'); %reversing y direction
     hold on; % Retain the contour plot
     plot(xGrid, yGrid, 'g+', 'MarkerSize', .5, 'LineWidth', .5); % green plus is the R_sub centre in C
     hold off; % Release the hold to prevent further additions
else
    % x_displacement is constant, 
    disp('x displacement all value are zero or constant no standard Deviation in data: Contour plot not generated.');
end
 
%ydisplacement plot

fprintf('The y_displacement array are:');
Y_Displacement= y_displacement

if std(y_displacement(:) ~= 0) 
     %see xdisplacement same line for understaing different line of this loop
     figure('Name','ydisplacement contour plot','NumberTitle','off');
     levels = linspace(min(y_displacement(:)), max(y_displacement(:)), 20); 
     contourf(xGrid, yGrid, y_displacement, levels, 'LineColor', 'r');
     colorbar; % Adds a color bar to indicate the displacement values
     title('Y Displacement Contour Plot');
     %colormap('hot'); % This can make small differences more visible colormap such as jet, parula, hot,
     xlabel('X-pixel');
     ylabel('Y -pixel');
     axis equal; % original size
     set(gca, 'YDir', 'reverse'); %reversing y direction
     hold on; % Retain the contour plot
     plot(xGrid, yGrid, 'g+', 'MarkerSize', .5, 'LineWidth', .5); % green plus is the R_sub centre in C
     hold off; % Release the hold to prevent further additions
else
    % y_displacement is constant, help eliminate error during validation
    disp('all x displacement point value is zero or constant or no Standard Deviation: Contour plot not generated.');
end

%final displacement plot

fprintf('The final_displacement array are:');
Final_Displacement= final_displacement

if std(final_displacement(:) ~= 0) 
     %  see xdisplacement for loop to  understand different line of this loop
     figure('Name','final_displacement contour plot','NumberTitle','off');
     levels = linspace(min(final_displacement(:)), max(final_displacement(:)), 10);  
     contourf(xGrid, yGrid, final_displacement, levels,'LineColor', 'r');
     colorbar; % Adds a color bar to indicate the displacement values
     title('final_displacement Contour Plot');
     %colormap('hot'); % This can make small differences more visible colormap such as jet, parula, hot,
     xlabel('X-pixel');
     ylabel('Y -pixel');
     axis equal; % original size
     set(gca, 'YDir', 'reverse'); %reversing y direction
     hold on; % Retain the contour plot
     plot(xGrid, yGrid, 'g+', 'MarkerSize', .5, 'LineWidth', .5); % green plus is the R_sub centre in C
     hold off; % Release the hold to prevent further additions
else
    % y_displacement is constant, 
    disp('final_displacement is constant  or no standart deviation : Contour plot not generated.');
end

%%%%%%%%%%%%%%%%%Step 6: Deformatino Gradient & jacobian %%%%%%%%%%%%

%Spatial Gradient of displacement field

[dx_dx, dx_dy]=gradient(x_displacement,grid_spaceing,grid_spaceing);%here grid_spaceing in x and y direction
[dy_dx, dy_dy]=gradient(y_displacement,grid_spaceing,grid_spaceing);

%initial 2*2 zeros matrix to hold deformation gradient(F)of each point

[a,b]=size(xGrid);
F=zeros(2,2,a,b);

for m=1:size(xGrid,1)
    for n=1:size(xGrid,2)
        F(:,:,a,b)= [dx_dx(i,j), dx_dy(i,j);
                 dy_dx(i,j), dy_dy(i,j)];
    end

end

fprintf('The Deformation Gradient array or Tensor (F) are:');
Tensor_F= F

%Jacobian Matrix

J=dx_dx.*dy_dy- dx_dy.*dy_dx;

fprintf('TheJacobian (J) are:');
Jacobian_J= J

%%%%%%%%%%%%%%step9: Lagrange strain & corresponding plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lagrange Strain(Exx,Eyy,Exy)

Exx=dx_dx+ 0.5*(dx_dx.^2+dy_dy.^2);
Eyy = dy_dy + 0.5*(dx_dy.^2 + dy_dy.^2);
Exy = 0.5*(dx_dy + dy_dx) + 0.5*(dx_dx.*dx_dy + dy_dx.*dy_dy);

%Exx plot

fprintf('The Exx strain array  are:');
Lagrange_Exx=Exx

if std(Exx(:) ~= 0) 
     %see xdisplacement comments to understand each line of this loop
     figure('Name','Exx contour plot','NumberTitle','off');
     levels = linspace(min(Exx(:)), max(Exx(:)), 5);
     contourf(xGrid, yGrid, Exx, levels, 'LineColor', 'r');
     colorbar; % Adds a color bar to indicate the displacement values
     title('Exx Contour Plot');
     %colormap('hot'); % This can make small differences more visible colormap such as jet, parula, hot,
     xlabel('X-pixel');
     ylabel('Y -pixel');
     axis equal; % original size
     set(gca, 'YDir', 'reverse'); %reversing y direction
     hold on; % Retain the contour plot
     plot(xGrid, yGrid, 'g+', 'MarkerSize', .5, 'LineWidth', .5); % green plus is the R_sub centre in C
     hold off; % Release the hold to prevent further additions
else
    % Exx is constant, 
    disp('Exx all value is constant or zeros no standard deviation: Contour plot not generated.');
end
 
%Eyy plot
fprintf('The Eyy strain array  are:');
Lagrange_Eyy=Eyy

if std(Eyy(:) ~= 0) 
     %see xdisplacement comments to understand each line of this loop
     figure('Name','Eyy contour plot','NumberTitle','off');
     levels = linspace(min(Eyy(:)), max(Eyy(:)), 5);
     contourf(xGrid, yGrid, Eyy, 'LineColor', 'r');
     colorbar; % Adds a color bar to indicate the displacement values
     title('Eyy Contour Plot');
     %colormap('hot'); % This can make small differences more visible colormap such as jet, parula, hot,
     xlabel('X-pixel');
     ylabel('Y -pixel');
     axis equal; % original size
     set(gca, 'YDir', 'reverse'); %reversing y direction
     hold on; % Retain the contour plot
     plot(xGrid, yGrid, 'g+', 'MarkerSize', .5, 'LineWidth', .5); % green plus is the R_sub centre in C
     hold off; % Release the hold to prevent further additions
else
    % Eyy  is constant, 
    disp('Eyy  all value is constant(zeros): Contour plot not generated.');
end

%Exy plot

fprintf('The Exy strain array  are:');
Lagrange_Exy=Exy

if std(Exy(:) ~= 0) 
     % if all data of the Exx is like of its first data 
     % i.e no strain - no graph.No error during valiadation process with same image, translation, rotation of the image
     figure('Name','Exy contour plot','NumberTitle','off');
     levels = linspace(min(Exy(:)), max(Exy(:)), 5);
     contourf(xGrid, yGrid, Exy, 'LineColor', 'r');
     colorbar; % Adds a color bar to indicate the displacement values
     title('Exy Contour Plot');
     %colormap('jet'); % This can make small differences more visible colormap such as jet, parula, hot,
     xlabel('X-pixel');
     ylabel('Y -pixel');
     axis equal; % original size
     set(gca, 'YDir', 'reverse'); %reversing y direction
     hold on; % Retain the contour plot
     plot(xGrid, yGrid, 'g+', 'MarkerSize', .5, 'LineWidth', .5); % green plus is the R_sub centre in C
     hold off; % Release the hold to prevent further additions
else
    % Exy  is constant, 
    disp('Exy  all value is constant(zeros) or  no standard deviation: Contour plot not generated.');
end

 