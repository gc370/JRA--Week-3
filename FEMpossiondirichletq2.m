function [approximationswithdirichlet,errormatrix] = FEMpossiondirichletq2(numberofpointsinx)

% This driver code is for solving the poisson equation; namely -u_xx - pi^2 u(x) = f(x)
                                                              % u(0) = 0,
                                                              % u(0.5) = 1
            
        % -------------- Inputs --------------- %
        
        %Boundary points coded in for this example.
        %xbeggining = 0 ;
        %xend = 1 ;
        %numberofpointsinx = 10; (your choice)

     actualsol = @(x) sin(pi*x);        

                
        %first set up the mesh
     xbeggining = 0;
     xend = 0.5;
     pointx = zeros(1,numberofpointsinx+1);
     dx = (xend - xbeggining)/(numberofpointsinx);
     diag = 2 - ((dx*pi)^2);

            for i=1:(numberofpointsinx+1)
                pointx(i) = xbeggining + (i-1)*dx;
            end

    % Now create the stiffness matrix
    
    matrixu = zeros(numberofpointsinx-1,numberofpointsinx-1);
    
    matrixu(1,1) = diag;
    matrixu(1,2) = -1;
    matrixu(numberofpointsinx-1,numberofpointsinx-2) = -1;
    matrixu(numberofpointsinx-1,numberofpointsinx-1) = diag;
    
             for i = 2:numberofpointsinx-2
                 matrixu(i,i-1) = -1;
                 matrixu(i,i) = diag;
                 matrixu(i,i+1) = -1;
            end
    
    % Now build the right hand side vector f.
    
    fknown = zeros(numberofpointsinx-1,1);
    
    fknown(numberofpointsinx-1) = 1;
    
    
    
    approximationsofu = inv(matrixu)*fknown; % Here matlab is inverting the matrix and multiplying it to fknown.
    approximationswithdirichlet = zeros(numberofpointsinx+1,1);
    approximationswithdirichlet(1) = 0;
    approximationswithdirichlet(numberofpointsinx+1) = 1;
    
            for i = 1:numberofpointsinx-1
                approximationswithdirichlet(i+1) =  approximationsofu(i);
            end
    
    correctsol = zeros(numberofpointsinx+1,1);
    for i = 1:numberofpointsinx+1
        
        correctsol(i) = actualsol(pointx(i));
        
    end
            
    errormatrix = abs(approximationswithdirichlet - correctsol);        
            
            
            
    %We now have our results, display as we wish
    
            for i = 1:numberofpointsinx+1
                disp("at x = " + pointx(i) + ", the value of U(" + pointx(i) + ") is " + approximationswithdirichlet(i))
            end
    
    
end