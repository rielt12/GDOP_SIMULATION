

x_1 = 0;
x_2 = 0;
x_3 = 0;


theta = [x_1; x_2; x_3];
y = [3/2; 1; -(10*pi-3)/3];

error =[
    
3*theta(1)-cos(theta(2)*theta(3))-y(1);
4*theta(1)^2-625*theta(2)^2+2*theta(2)-y(2);
exp(-theta(1)*theta(2))+20*theta(3)-y(3);

]; 


J = [3, sin(theta(2)*theta(3))*theta(3),sin(theta(2)*theta(3))*theta(2);
    8*theta(1), -1250*theta(2)+2, 0 ;
    -theta(2)*exp(-theta(1)*theta(2)),-theta(1)*exp(-theta(1)*theta(2)),20];


    
cost = 0.5*norm(error)^2;

iter = 1;
while(cost> 1e-9)
    cost = 0.5*norm(error)^2
    Cost(iter,1) = cost;

    delta_theta = J'*error;

    gamma = line_search(theta,delta_theta, 10000, y);
    Gamma(iter,1) = gamma;

    theta =  theta- gamma*delta_theta; 


J = [3, sin(theta(2)*theta(3))*theta(3),sin(theta(2)*theta(3))*theta(2);
    8*theta(1), -1250*theta(2)+2, 0 ;
    -theta(2)*exp(-theta(1)*theta(2)),-theta(1)*exp(-theta(1)*theta(2)),20];

    error =[
    
3*theta(1)-cos(theta(2)*theta(3))-y(1);
4*theta(1)^2-625*theta(2)^2+2*theta(2)-y(2);
exp(-theta(1)*theta(2))+20*theta(3)-y(3);

]; 


    iter = iter + 1;
end

figure(1)
plot(Gamma)


figure(2)
plot(Cost)

