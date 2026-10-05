function values = physical_rate_design_values(type,x)
switch upper(type)
    case '4L'
        values = [10.^x(1),10.^x(2),10.^x(3),x(4)];
    case '5L'
        values = [10.^x(1),10.^x(2),10.^x(3),10.^x(4),x(5)];
    otherwise
        error('Unknown receiver type: %s',type);
end
end
