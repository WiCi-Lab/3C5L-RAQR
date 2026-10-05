function x = linear_level_crossing(x1,x2,y1,y2,target)
if abs(y2-y1)<=realmin
    x=(x1+x2)/2;
else
    x=x1+(target-y1)*(x2-x1)/(y2-y1);
end
end
