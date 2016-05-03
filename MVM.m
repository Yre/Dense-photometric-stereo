clear global;
global vertice;


icosahedron();

%% Uniform resampling
datapath = './data/data04/';
lightvec = load([datapath 'lightvec.txt']);
[IDX,~]=knnsearch(vertice,lightvec);
[new_vertice,~,reverse_idx] = unique(IDX);

num_direction = size(new_vertice,1);

image_files = dir( fullfile( datapath, '*.bmp') ); 
num_images = length(image_files);
I=imread(fullfile(datapath, image_files(1).name));

normal_images = double(zeros([num_direction size(I)]));
weight = zeros(num_direction,1);

for i=1:num_images
    file=fullfile(datapath, image_files(i).name);
    Img = imread(file);
    normal_images(reverse_idx(i),:,:,:) = normal_images(reverse_idx(i),:,:,:) + double([0 (lightvec(i,:)*vertice(IDX(i),:)')*Img]);
    weight(reverse_idx(i)) = weight(reverse_idx(i)) + dot(lightvec(i),vertice(IDX(i)));
end

for i=1:num_direction
    normal_images(i,:,:,:) = normal_images(i,:,:,:)/weight(i);
end

%% Choose Denominator

for i=1:num_direction
    grey_images(i,:,:) = rgb2grey(normal_images(i,:,:,:));
end








