function [isDisallowedBoundary,x_DC,y_DC] = DCClustersToGrid(NumDCs,R_DC,numberOfClusters,domainBoundary,disallowedBoundary)
    
    
    isDisallowedBoundary = zeros(size(disallowedBoundary.X,1),size(disallowedBoundary.Y,2));
    
    [x_DC,y_DC] = DC_Clusters(NumDCs,R_DC+1.0,numberOfClusters,domainBoundary);

    if false
        figure;
        scatter(x_DC,y_DC,'filled')
        axis([domainBoundary.x_min domainBoundary.x_max domainBoundary.y_min domainBoundary.y_max])
        drawnow;
    end
    
    for jj=1:size(disallowedBoundary.Y,1)
        for ii=1:size(disallowedBoundary.X,2)
            for kk=1:length(x_DC)
                if ( (ii-x_DC(kk))^2 + (jj-y_DC(kk))^2 <= R_DC^2 )
                    isDisallowedBoundary(jj,ii) = 1;
                end
            end
           
        end
    end
end
