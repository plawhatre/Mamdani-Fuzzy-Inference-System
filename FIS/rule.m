%Rule base
function [ r3 ] = rule( r1,r2 )
if r1==1 && r2==1
    r3=1;
end
if r1==1 && r2==2
    r3=1;
end
if r1==1 && r2==3
    r3=2;
end
if r1==2 && r2==1
    r3=2;
end
if r1==2 && r2==2
    r3=2;
end
if r1==2 && r2==3
    r3=3;
end
if r1==3 && r2==1
    r3=3;
end
if r1==3 && r2==2
    r3=3;
end
if r1==3 && r2==3
    r3=4;
end
if r1==4 && r2==1
    r3=4;
end
if r1==4 && r2==2
    r3=4;
end
if r1==4 && r2==3
    r3=4;
end
end

