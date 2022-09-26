function phi_new = RK2SecondOBackward(phi_old,t_end,Delta_t,CFL)

N_points=size(phi_old,2);
t = 0;
phi_new=phi_old;


while (t < t_end)
    
        phi_old = phi_new;
     
        phi_intermediate(1) = phi_old(1) - 0.5 * CFL * (1.5*phi_old(1) - ...
            2*phi_old(N_points-1) + 0.5*phi_old(N_points-2));
        
        phi_intermediate(2) = phi_old(2)-0.5* CFL * (1.5*phi_old(2) - ...
            2*phi_old(1) + 0.5*phi_old(N_points-1));
        
        for i=3:N_points
            phi_intermediate(i) = phi_old(i)-0.5* CFL * (1.5*phi_old(i) - ...
                2*phi_old(i-1) + 0.5*phi_old(i-2));
        end

        
        phi_new(1) = phi_old(1)- CFL * (1.5*phi_intermediate(1) - ...
            2*phi_intermediate(N_points-1) + 0.5*phi_intermediate(N_points-2));
        
        phi_new(2) = phi_old(2)- CFL * (1.5*phi_intermediate(2) - ...
            2*phi_intermediate(1) + 0.5*phi_intermediate(N_points-1));
        
        for i=3:N_points
            phi_new(i) = phi_old(i)-CFL * (1.5*phi_intermediate(i) - ...
                2*phi_intermediate(i-1) + 0.5*phi_intermediate(i-2));
        end
    
        
    t=t+Delta_t;
    
end


end