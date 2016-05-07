clear global;
<<<<<<< HEAD
global vertice;


icosahedron();

%% Uniform resampling
datapath = './data/data04/';
=======
addpath(genpath('./lib/'));
addpath(genpath('./util/'));

vertice = icosahedron(4);
vertice = vertice(vertice(:,3)>=0,:);

% fclose(fileID);
%% Uniform resampling
datapath = '../data/data02/';
>>>>>>> 2b7efc3b742c2d428ca87e21aa31f800e1888b3d
lightvec = load([datapath 'lightvec.txt']);
[IDX,~]=knnsearch(vertice,lightvec);
[new_vertice,~,reverse_idx] = unique(IDX);

num_direction = size(new_vertice,1);

image_files = dir( fullfile( datapath, '*.bmp') ); 
num_images = length(image_files);
I=imread(fullfile(datapath, image_files(1).name));

normal_images = double(zeros([size(I) num_direction]));
weight = zeros(num_direction,1);

for i=1:num_images
    file=fullfile(datapath, image_files(i).name);
<<<<<<< HEAD
    Img = imread(file);
    normal_images(reverse_idx(i),:,:,:) = normal_images(reverse_idx(i),:,:,:) + double([0 (lightvec(i,:)*vertice(IDX(i),:)')*Img]);
    weight(reverse_idx(i)) = weight(reverse_idx(i)) + dot(lightvec(i),vertice(IDX(i)));
=======
    tic_toc_print('Interpolate images %d / %d\n', i, num_images);
    Img = double(imread(file));
    i_w = lightvec(i,:)*vertice(IDX(i),:)';
    normal_images(:,:,:,reverse_idx(i)) = normal_images(:,:,:,reverse_idx(i)) + Img*i_w;
    weight(reverse_idx(i)) = weight(reverse_idx(i)) + i_w;
>>>>>>> 2b7efc3b742c2d428ca87e21aa31f800e1888b3d
end

for i=1:num_direction
    normal_images(:,:,:,i) = normal_images(:,:,:,i)/weight(i);
end
lightvec = vertice(new_vertice,:);
%% Choose Denominator
<<<<<<< HEAD

for i=1:num_direction
    grey_images(i,:,:) = rgb2grey(normal_images(i,:,:,:));
=======
grey_images=double(zeros(size(I,1),size(I,2),num_direction));
for i=1:num_direction
%     grey_images(i,:,:) = rgb2gray(squeeze(normal_images(:,:,:,i))/255);
  grey_images(:,:,i) = 0.2989 * normal_images(:,:,1,i) + ...
    0.5870 * normal_images(:,:,2,i) + 0.1140 * normal_images(:,:,3,i);

>>>>>>> 2b7efc3b742c2d428ca87e21aa31f800e1888b3d
end
pix_Rank = zeros(size(grey_images));
AA = zeros(num_direction,1);

<<<<<<< HEAD




=======
for i=1:size(I,1)
    for j=1:size(I,2)
        pixel=squeeze(grey_images(i,j,:));
        [~,id_]=sort(pixel);
        AA(id_)=1:num_direction;
        pix_Rank(i,j,:)=AA;
    end
end
L = 0.7 * num_direction;
H = 0.9 * num_direction;
k_L = zeros(num_direction,1);
r_L = zeros(num_direction,1);
for i=1:num_direction
    threshold = pix_Rank(:,:,i)>L;
    k_L(i) = sum(sum(threshold(:)));
    r_L(i) = mean2(grey_images(:,:,i).*threshold);
   
end
k_L =  k_L .*(r_L<H);
[~,denominator_idx]=max(k_L);
disp(denominator_idx);

%% Local normal estimation 
% ratio_images = grey_images ./ repmat(grey_images(denominator_idx,:,:),num_direction,1);
ratio_images = double(grey_images);
%if the denominator not removed?
init_norm = zeros(size(I,1),size(I,2),3);
light_deno = lightvec(denominator_idx,:);
light_ratio = lightvec(1:denominator_idx-1,:);
light_ratio = [light_ratio; lightvec(denominator_idx+1:end,:)];
fileID=fopen('init_norm.txt','w');
for i=1:size(I,1)
    for j=1:size(I,2)
        I_deno = ratio_images(i,j,denominator_idx);
        I_ratio = [squeeze(ratio_images(i,j,1:denominator_idx-1)); squeeze(ratio_images(i,j,denominator_idx+1:end))];  
        A = light_ratio .* I_deno  - I_ratio * light_deno; 
        [~,~,N] = svd(A,0);
        if (N(3,3)>0)
            init_norm(i,j,:) = N(:,3);
        else
            init_norm(i,j,:) = -N(:,3);
        end
        fprintf(fileID,'%.3f %.3f %.3f\n',init_norm(i,j,:));
    end
end
>>>>>>> 2b7efc3b742c2d428ca87e21aa31f800e1888b3d

fclose(fileID);
figure, imshow(init_norm);

%% Refine
% opt_matrix = graphcut_refine(init_norm);

sigma = 0.6;
lambda = 0.5;

vertice = icosahedron(5);
vertice = vertice(vertice(:, 3) > 0,:);
norm_vec = reshape(init_norm, [],3);
[IDX, ~] = knnsearch(vertice, norm_vec);
[new_vertice, ~, ~] = unique(IDX);
vertice = vertice(new_vertice, :);

E_data = pdist2(vertice, norm_vec);
E_smoothness = lambda * log(1 + pdist2(vertice, vertice) / (2 * sigma * sigma));
%convert data type in energy as integer by divide smallest non-zero number
 % ???
E_data = int32(E_data * 10000);
E_smoothness = int32(E_smoothness * 10000);


edge_num = (size(I,1)-1)*size(I,2) + (size(I,1))*(size(I,2)-1);
Si=zeros(edge_num,1);
Sj=zeros(edge_num,1);
Sv=ones(edge_num,1);
cnt = 0;
%vertical
for i=1:size(I,1)-1
    for j=1:size(I,2)
        cnt=cnt+1;
        Si(cnt)=(i-1)*size(I,2)+j;
        Sj(cnt)=i*size(I,2)+j;
    end
end
%horizontal
for i=1:size(I,1)
    for j=1:size(I,2)-1
        cnt=cnt+1;
        Si(cnt)=(i-1)*size(I,2)+j;
        Sj(cnt)=(i-1)*size(I,2)+j+1;
    end
end

s = size(init_norm);
labels = vertice;

L = size(labels,1); 

neighbor_matrix = sparse(Si,Sj,Sv,s(1)*s(2),s(1)*s(2));

h=GCO_Create(s(1)*s(2),L);
GCO_SetDataCost(h, E_data);
GCO_SetSmoothCost(h, E_smoothness);
GCO_SetNeighbors(h, neighbor_matrix);
GCO_Expansion(h);
opt_label = GCO_GetLabeling(h);
GCO_Delete(h);


normsOPT1D = norm_vec(opt_label,:);

opt_norm = zeros(size(init_norm));
for j = 1:size(I,2)
  for i = 1:size(I,1)
    opt_norm(i,j,:) = normsOPT1D((j-1)*size(I,1)+i, :);
  end
end
opt_norm = vertice(opt_label, :);
opt_norm = reshape(opt_norm, [size(I,1), size(I,2), 3]);

normsOPT = opt_norm;
figure, imshow((-1/sqrt(3) * normsOPT(:,:,1) + 1/sqrt(3) * normsOPT(:,:,2) + 1/sqrt(3) * normsOPT(:,:,3)) / 1.1);


%% reconstruct

opt_norm = init_norm;

map_height = size(opt_norm, 1);
map_width = size(opt_norm, 2);
slant = zeros(map_height, map_width);
tilt = zeros(map_height, map_width);

for i = 1:map_height
  for j = 1:map_width
    % note need to figure out why we need different coordinate system
    vec = squeeze(opt_norm(map_height+1-i, j, :));
    x = vec(1); y = vec(2); z = vec(3);
    dzdx = -x / z; dzdy = -y / z;
    [slant(i,j), tilt(i,j)] = grad2slanttilt(dzdx,dzdy);
  end
end

depth_map = shapeletsurf(slant, tilt, 6, 3, 2);
figure, surf(depth_map);
