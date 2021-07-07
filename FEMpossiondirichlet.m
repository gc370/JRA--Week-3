function [approximationsofu] = FEMpossiondirichlet(f,g1,g2,numberofpointsinx,xbeggining,xend)

% This driver code is for solving the poisson equation; namely -u_xx = f(x)
            
        % -------------- Inputs --------------- %
        
        %g1= 1; -> first dirichlet boundary condition
        %g2= exp(-1); -> Second dirichlet boundary condition
        %xbeggining = 0 ;
        %xend = 1 ;
        %f= @(x) -exp(-x);
        %numberofpointsinx = 10; (your choice)

        

                
        %first set up the mesh


     pointx = zeros(1,numberofpointsinx+1);
     dx = (xend - xbeggining)/(numberofpointsinx);

            for i=1:(numberofpointsinx+1)
                pointx(i) = xbeggining + (i-1)*dx;
            end

    % Now create the stiffness matrix
    
    matrixu = zeros(numberofpointsinx-1,numberofpointsinx-1);
    
    matrixu(1,1) = 2;
    matrixu(1,2) = -1;
    matrixu(numberofpointsinx-1,numberofpointsinx-2) = -1;
    matrixu(numberofpointsinx-1,numberofpointsinx-1) = 2;
    
             for i = 2:numberofpointsinx-2
                 matrixu(i,i-1) = -1;
                 matrixu(i,i) = 2;
                 matrixu(i,i+1) = -1;
            end
    
    % Now build the right hand side vector f.
    
    fknown = zeros(numberofpointsinx-1,1);
    
            for i=1:numberofpointsinx-1
                fknown(i) = (dx)^2 * f(pointx(i+1));      
            end
            
    fknown(1) = fknown(1) + g1;                                      
    fknown(numberofpointsinx-1) = fknown(numberofpointsinx-1) + g2;
    
    
    
    approximationsofu = inv(matrixu)*fknown; % Here matlab is inverting the matrix and multiplying it to fknown.
    approximationswithdirichlet = zeros(numberofpointsinx+1,1);
    approximationswithdirichlet(1) = g1;
    approximationswithdirichlet(numberofpointsinx+1) = g2;
    
            for i = 1:numberofpointsinx-1
                approximationswithdirichlet(i+1) =  approximationsofu(i);
            end
    
    %We now have our results, display as we wish
    
            for i = 1:numberofpointsinx+1
                disp("at x = " + pointx(i) + ", the value of U(" + pointx(i) + ") is " + approximationswithdirichlet(i))
            end
    
    
end

