clear all
MD = importdata('neighbor_vector_L20_np20.dat');
%num_particles = 10;
list_length = length(MD.data(:,1)) - 1;

particles = unique(MD.data(1:end-1,1),'stable');
num_particles = length(particles);
num_frames = list_length / num_particles;
k = 1;
% 上面是得出总共有几个particle（num_particles)
%和循环了几次/要画几张图（num_frames)


% 把每个frame所需要的x，y坐标集起来
for i = 1:num_frames
    %从第一个frame开始
    l = 1;
  %  while (k <= list_length)
        for j = k:(k+num_particles)
        
            tempx(l) = MD.data(j,2);
            tempy(l) = MD.data(j,3);
            l = l+1;
       
        end
   % end
    frames(1,:,i) = tempx;
    frames(2,:,i) = tempy;
%     frames(1,:,i) = tempx;
%     frames(2,:,i) = [frames(i).y, MD.data(k,3)];
    k = k + num_particles;
    l = 0;
end

xmax = max(MD.data(:,2)) + 1;
ymax = max(MD.data(:,3)) + 1;
xmin = min(MD.data(:,2)) - 1;
ymin = min(MD.data(:,3)) - 1;
[x,y,z] = size(frames);
%num_frames = len / num_particles;

for jello = 1:num_frames
    for point = 1:y
        plot(frames(1,point,jello),frames(2,point,jello),'.','MarkerSize',20)
        hold on
    end
    grid on
    xlim([xmin,xmax]);
    ylim([ymin,ymax]);
    title('Particles Movement','FontSize',20);
    hold off
    F(jello) = getframe; 
end

writerObj.FrameRate = 0.5;
writerObj = VideoWriter('neighbor_vector_L20_np20.avi');

open(writerObj);
writeVideo(writerObj, F)
close(writerObj);

