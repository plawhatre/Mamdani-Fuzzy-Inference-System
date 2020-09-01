function [ cx,cy ] = cntr( x,y )
a=0;
for i=1:(length(x)-1)
    a=a+(1/2)*(x(i)*y(i+1)-x(i+1)*y(i));
end
cx=0;
cy=0;
for i=1:(length(x)-1)
    cx=cx+(x(i)+x(i+1))*(x(i)*y(i+1)-x(i+1)*y(i))/(6*a);
    cy=cy+(y(i)+y(i+1))*(x(i)*y(i+1)-x(i+1)*y(i))/(6*a);
end

end

