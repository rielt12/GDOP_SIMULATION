function Kost = compute_cost(theta , y)
error =[
    
3*theta(1)-cos(theta(2)*theta(3))-y(1);
4*theta(1)^2-625*theta(2)^2+2*theta(2)-y(2);
exp(-theta(1)*theta(2))+20*theta(3)-y(3);

]; 


Kost = 0.5*norm(error)^2;

end