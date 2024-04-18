function [s, s_dot, s_dotdot]=trap_profile(t_i ,t_f, f_s, p_i, p_f)
    arguments % default argument
        t_i = 0.0
        t_f = 0.0
        f_s = 1000
        p_i = 0.0
        p_f = 1.0
    end

s_i = 0;
s_f = norm(p_f-p_i);

delta_t = t_f-t_i;

s_c_dot = 1.5*abs(s_f-s_i) / delta_t;
tc = (s_i-s_f+s_c_dot*delta_t) / s_c_dot;
s_c_dotdot = s_c_dot/tc;

numero_campioni = f_s*delta_t;


t = linspace(0, delta_t, numero_campioni);
s = zeros(1, numero_campioni);
s_dot = zeros(1, numero_campioni);
s_dotdot = zeros(1, numero_campioni);

for k = 1:numero_campioni
    if t(k) <= tc
        s(k) = s_i+(1/2)*s_c_dotdot*t(k).^2;
        s_dot(k) = s_c_dotdot*t(k);
        s_dotdot(k) = s_c_dotdot;
    elseif ( t(k)>tc && t(k) <= delta_t-tc)
            s(k) = s_i+s_c_dotdot*tc*(t(k)-tc/2);
            s_dot(k) = s_c_dotdot*tc;
            s_dotdot(k) = 0;
    elseif ( t(k) > delta_t-tc && t(k) <= delta_t)
            s(k)=s_f-(1/2)*s_c_dotdot*(delta_t-t(k)).^2;
            s_dot(k)=s_c_dotdot*(delta_t-t(k));
            s_dotdot(k)=-s_c_dotdot;
    end
end

end