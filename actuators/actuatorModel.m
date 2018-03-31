function [ Tm, Trw, wrw_new ] = actuatorModel( Tctrl, B,  wrw )

mu = cross(Tctrl, B);
musat = 0.015;
if max(mu) <= musat
    mu_m = mu;
else
    mu_unit = mu/norm(mu);
    mu_max = max(mu);T
    mu_m = musat/mu_max*mu_unit;
end
Tm = cross(mu_m, B);

Trem = Tctrl - Tm;
if Trem == zeros(3,1)
    
Irw = 3.534e-3/(7500*2*pi/60);
dwrw = Trem/Irw; 
%somehow cap the increase in reaction wheel speed

wrw_new = wrw+dwrw;
Trw = Irw*dwrw;

end

end