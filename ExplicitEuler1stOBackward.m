function phi_new = ExplicitEuler1stOBackward(phi_old,t_end,Delta_t,CFL)

N_points=size(phi_old,2);
t = 0;
phi_new=phi_old;


while (t < t_end)
    
    phi_old = phi_new;
    
    phi_new(1) = phi_old(1)- CFL * (phi_old(1) - phi_old(N_points-1));
    for i=2:N_points
        phi_new(i) = phi_old(i)- CFL * (phi_old(i) - phi_old(i-1));
    end
    
    t=t+Delta_t;
    
end


end