function y = cosfit_rad(params,r)

lag = params(1,1);
int = params(1,2);
amp = params(1,3);

y = int + amp.*(cos(r + lag));

end