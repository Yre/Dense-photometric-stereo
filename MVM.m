clear global;


vertice = icosahedron(4);
fileID = fopen('icosahedron_vertices.txt','w');
vertice = vertice(vertice(:,3)>=0,:);

fprintf(fileID,'%.4f %.4f %.4f\n',vertice);
fclose(fileID);
%% Uniform resampling
datapath = '../data/data04/';
lightvec = load([datapath 'lightvec.txt']);
[IDX,~]=knnsearch(vertice,lightvec);
[new_vertice,~,reverse_idx] = unique(IDX);

num_direction = size(new_vertice,1);

image_files = dir( fullfile( datapath, '*.bmp') ); 
num_images = length(image_files);
I=imread(fullfile(datapath, image_files(1).name));

normal_images = double(zeros([num_direction size(I)]));
weight = zeros(num_direction,1);
tmp_img=zeros([1 size(I)]);

for i=1:num_images
    file=fullfile(datapath, image_files(i).name);
    Img = imread(file);
    tmp_img(1,:,:,:) = (Img .* dot(lightvec(i,:),vertice(IDX(i),:)));
    normal_images(reverse_idx(i),:,:,:) = normal_images(reverse_idx(i),:,:,:) + tmp_img;
    weight(reverse_idx(i)) = weight(reverse_idx(i)) + dot(lightvec(i),vertice(IDX(i)));
end
lightvec = vertice(new_vertice,:);
for i=1:num_direction
    normal_images(i,:,:,:) = normal_images(i,:,:,:)/weight(i);
end

%% Choose Denominator
grey_images=double(zeros(num_direction,size(I,1),size(I,2)));
for i=1:num_direction
    grey_images(i,:,:) = rgb2gray(squeeze(normal_images(i,:,:,:) / 255));
end

pixel=zeros(size(I,1),size(I,2));
for i=1:size(I,1)
    for j=1:size(I,2)
        pixel=squeeze(grey_images(:,i,j));
        [~,id_]=sort(pixel);
        AA(id_)=1:num_direction;
        pix_Rank(i,j,:)=AA;
    end
end

L = 0.7 * num_direction;
H = 0.9 * num_direction;

for i=1:num_direction
    threshold = pix_Rank(:,:,i)>L;
    k_L(i) = sum(sum(threshold));
    r_L(i) = mean2(pix_Rank(:,:,i).*threshold);
end

k_L = (r_L<H) .* k_L;
[~,denominator_idx]=max(k_L);
disp(denominator_idx);

%% Local normal estimation 
ratio_images = grey_images ./ repmat(grey_images(denominator_idx,:,:),num_direction,1);

%if the denominator not removed?
init_norm = zeros(size(I,1),size(I,2),3);
light_deno = lightvec(denominator_idx,:);
light_ratio = lightvec(1:denominator_idx-1,:);
light_ratio = [light_ratio; lightvec(denominator_idx+1,:)];

for i=1:size(I,1)
    for j=1:size(I,2)
        I_deno = ratio_images(denominator_idx,i,j);
        I_ratio = squeeze( [ratio_images(1:denominator_idx-1,i,j); ratio_images(denominator_idx+1,i,j)]);  
        A = light_ratio .* I_deno  - I_ratio * light_deno; 
        [~,~,N] = svd(A,0);
        if (N(3,3)>0)
            init_norm(i,j,:) = N(:,3);
        else
            init_norm(i,j,:) = -N(:,3);
        end
    end
end



