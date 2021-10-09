% points in the world coordinate 
P = [ 4 0 4;
      4 0 6;
      3 0 7;
      6 0 8;
      2 0 2;
      6 0 2;
      0 1 9;
      0 3 8;
      0 4 5;
      0 6 7;
      0 2 2;
      0 7 3];
      
% points in 2D         
p = [ 219 269;
      213 198;
      244 160;
      133 120;
      280 321;
      150 354;
      351 88;
      402 119;
      432 233;
      505 155;
      377 323;
      538 340];
      
n=12;

% normalisation of the data
% Real world
P_center = mean(P);      
P0 = P-P_center;
avg_dist_P0 = 0;
for i=1:n
  avg_dist_P0 = avg_dist_P0 + norm(P0(i,:),'fro');
end
avg_dist_P0 = avg_dist_P0/n;
scale1 = sqrt(3)/avg_dist_P0;
P_norm = scale1*P0;
P_norm(:,4) = ones(n,1);
P_scale = diag([scale1 scale1 scale1 1]);
P_shift = eye(4);
P_shift(:,4)=[-P_center 1];
P_trans = P_scale*P_shift;

% Image plane
p_center = mean(p);
p0 = p - p_center;
avg_dist_p0 = 0;
for i=1:n
    avg_dist_p0 = avg_dist_p0 + norm(p0(i,:),'fro');
end
avg_dist_p0 = avg_dist_p0/n;
scale2 = sqrt(2)/avg_dist_p0;
p_norm = scale2*p0;
p_norm(:,3) = ones(n,1);
p_scale = diag([scale2 scale2 1]);
p_shift = eye(3);
p_shift(:,3)=[-p_center 1];
p_trans = p_scale*p_shift;


p(:,3) = ones(n,1);
P(:,4) = ones(n,1);

% DLT method
O = [ 0 0 0 0 ];
Q = zeros(2*n,12);
for i=1:n
  Q(2*i-1,:) = [ P_norm(i,:) O -p_norm(i,1)*P_norm(i,:) ];
  Q(2*i,:) = [ O P_norm(i,:) -p_norm(i,2)*P_norm(i,:) ];
end
[V,D] = eig(Q.'*Q);
M_norm = reshape(V(:, 1), [], 3)';
M = (inv(p_trans)*M_norm)*P_trans;
Project_Pts = ((M*P')');
Project_Pts = Project_Pts./Project_Pts(:, 3);

% error in projected points and actual 2D points
E = Project_Pts-p;
SQE  = E.^2;
MSE  = mean(SQE(:));
RMSE = sqrt(MSE);

% Plotting Actual points and projected points
plot(p(:,1),p(:,2),'ro');
hold on
plot(Project_Pts(:,1),Project_Pts(:,2),'+');
legend('Actual Points','Projected Points');
hold off


% calcuating intrinsic and extrinsic parameters
A = M(:,1:3);
R = zeros(3,3);
K = zeros(3,3);
rho = -1/norm(A(3,:),2);
R(3,:) = rho*A(3,:);
X0 = rho*rho*(A(1,:)*A(3,:).')
Y0 = rho*rho*(A(2,:)*A(3,:).')
cross1 = cross(A(1,:),A(3,:));
cross2 = cross(A(2,:),A(3,:));
n_cross1 = norm(cross1,2);
n_cross2 = norm(cross2,2);
theta = acos(-cross1/n_cross1*cross2.'/n_cross2)
alpha = rho*rho*n_cross1*sin(theta)
Beta = rho*rho*n_cross2*sin(theta)
R(1,:) = cross2/n_cross2;
R(2,:) = cross(R(3,:),R(1,:))
K(1,:) = [alpha -alpha*cot(theta) X0];
K(2,:) = [0 Beta/sin(theta) Y0];
K(3,:) = [0 0 1]
t = (rho*inv(K)*M(:,4)')'
x0 = -inv(R)*t
