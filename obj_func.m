function sqloss=obj_func(w,theta,y)
theta_w = theta*w;
sqloss = mean((theta_w-y)'*(theta_w-y));
end