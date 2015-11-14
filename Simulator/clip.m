hold on
L1 = data.ship(1,1).Lpp;
B1 = data.ship(1,1).B;
L2 = data.ship(1,2).Lpp;
B2 = data.ship(1,2).B;
for kp1 = 1:10:lt
    [l1(1),l1(2),l1(3),l1(4),l1(5)] = shipdraw(variable.ship(1,1).eta(1,kp1),variable.ship(1,1).eta(2,kp1),variable.ship(1,1).eta(6,kp1),L1,B1,1,-1,[0 0 1]);
    [l2(1),l2(2),l2(3),l2(4),l2(5)] = shipdraw(variable.ship(1,2).eta(1,kp1),variable.ship(1,2).eta(2,kp1),variable.ship(1,2).eta(6,kp1),L2,B2,1,-1,[1 0 0]);
    str_t = ['t = ',num2str(tsim(kp1)),'s'];
    d = sqrt((variable.ship(1,1).eta(1,kp1)-variable.ship(1,2).eta(1,kp1))^2+(variable.ship(1,1).eta(2,kp1)-variable.ship(1,2).eta(2,kp1))^2);
    str_d = ['d = ',num2str(d),'m'];
    t = text(variable.ship(1,1).eta(2,kp1)-100,variable.ship(1,1).eta(1,kp1),str_t);
    dist = text(variable.ship(1,1).eta(2,kp1)-100,variable.ship(1,1).eta(1,kp1)-20,str_d);
    pause(0.1);
    for kp2 = 1:5
        set(l1(kp2),'visible','off')
        set(l2(kp2),'visible','off')
        set(t,'visible','off')
        set(dist,'visible','off')
    end
end


