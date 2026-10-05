function p2 = decode_rate_design(type,p,x)
p2 = p;

switch upper(type)
    case '4L'
        if numel(x)~=4, error('4L rate design requires four variables.'); end
        p2.E_LO_RF_4L = 10.^x(1);
        p2.P_p_4L = 10.^x(2)*1e-3;
        p2.P_c_4L = 10.^x(3)*1e-3;
        p2.detuning_RF_4L = x(4)*2*pi*1e6;

    case '5L'
        if numel(x)~=5, error('5L rate design requires five variables.'); end
        p2.E_LO_RF_5L = 10.^x(1);
        p2.P_p_5L = 10.^x(2)*1e-3;
        p2.P_d_5L = 10.^x(3)*1e-3;
        p2.P_c_5L = 10.^x(4)*1e-3;
        p2.detuning_RF_5L = x(5)*2*pi*1e6;

    otherwise
        error('Unknown receiver type: %s',type);
end
end
