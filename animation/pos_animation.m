MD = importdata('test.dat');
num_particles = 10;
len = length(MD.data(:,1));

a = 1;
b = 1;
c = 1;
d = 1;
e = 1;
f = 1;
g = 1;
h = 1;
ai = 1;
j = 1;

% GET POSITIONS FOR EACH PARTICLE
for i = 1:len
    
    if (MD.data(i,1) == 0)
        particle0.x(a) = MD.data(i,2);
        particle0.y(a) = MD.data(i,3);
        a = a + 1;
    elseif (MD.data(i,1) == 1)
        particle1.x(b) = MD.data(i,2);
        particle1.y(b) = MD.data(i,3);
        b = b + 1;
    elseif (MD.data(i,1) == 2)
        particle2.x(c) = MD.data(i,2);
        particle2.y(c) = MD.data(i,3);
        c = c + 1;
    elseif (MD.data(i,1) == 3)
        particle3.x(d) = MD.data(i,2);
        particle3.y(d) = MD.data(i,3);
        d = d + 1;
    elseif (MD.data(i,1) == 4)
        particle4.x(e) = MD.data(i,2);
        particle4.y(e) = MD.data(i,3);
        e = e + 1;
    elseif (MD.data(i,1) == 5)
        particle5.x(f) = MD.data(i,2);
        particle5.y(f) = MD.data(i,3);
        f = f + 1;
    elseif (MD.data(i,1) == 6)
        particle6.x(g) = MD.data(i,2);
        particle6.y(g) = MD.data(i,3);
        g = g + 1;
    elseif (MD.data(i,1) == 7)
        particle7.x(h) = MD.data(i,2);
        particle7.y(h) = MD.data(i,3);
        h = h + 1;
    elseif (MD.data(i,1) == 8)
        particle8.x(ai) = MD.data(i,2);
        particle8.y(ai) = MD.data(i,3);
        ai = ai + 1;
    elseif (MD.data(i,1) == 9)
        particle9.x(j) = MD.data(i,2);
        particle9.y(j) = MD.data(i,3);
        j = j + 1;
    end        
end

num_frames = len / num_particles;

for jello = 1:num_frames
    plot(particle0.x(jello),particle0.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle1.x(jello),particle1.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle2.x(jello),particle2.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle3.x(jello),particle3.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle4.x(jello),particle4.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle5.x(jello),particle5.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle6.x(jello),particle6.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle7.x(jello),particle7.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle8.x(jello),particle8.y(jello),'.','MarkerSize',20)
    hold on
    plot(particle9.x(jello),particle9.y(jello),'.','MarkerSize',20)
    xlim([-5,12]);
    ylim([-5,12]);
    title('particle movement');
    F(jello) = getframe; 
end

writerObj = VideoWriter('particle_movement.avi');
open(writerObj);
writeVideo(writerObj, F)
close(writerObj);




