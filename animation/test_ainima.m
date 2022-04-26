time = 10;
for t = 1:time
   fplot(@(x) sin(x*50/t),[0,2*pi]);  % plot
   ylim([-1,1]);                      % guarantee consistent height
   F(t) = getframe;                   % capture it
end

writerObj = VideoWriter('particle_movement.avi');
open(writerObj);
writeVideo(writerObj, F)
close(writerObj);