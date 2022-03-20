function [GrPRS, GsPRS, GrRad, GsRad, R] = calc_slice(GrPRS, theta)
    R0 = [0,1,0;
         -1,0,0;
          0,0,1];
    iR0 =[0,-1,0;
          1, 0,0;
          0, 0,1];
    Gr0 = R0(2,:);
    Gs0 = R0(3,:);
    GsPRS = zeros(size(GrPRS));
    GpPRS = zeros(size(GrPRS));
    
    GrRad = zeros(size(GrPRS));
    GsRad = zeros(size(GrPRS));
    GpRad = zeros(size(GrPRS));
    R = zeros(size(GrPRS,1), 3, 3);
    for ii = 1:size(GrPRS,1) 
        if Gs0*GrPRS(ii,:)' > 0
            GsPRS(ii,:) = cross(Gr0, GrPRS(ii,:));
        else
            GsPRS(ii,:) = cross(Gs0, GrPRS(ii,:));
        end
        
        GsPRS(ii,:) = GsPRS(ii,:) / norm(GsPRS(ii,:));
        GpPRS(ii,:) = cross(GrPRS(ii,:), GsPRS(ii,:));
        GpPRS(ii,:) = GpPRS(ii,:) / norm(GpPRS(ii,:));
        
        GrRad(ii,:) = iR0 * reshape(GrPRS(ii,:), [3,1]);
        GsRad(ii,:) = iR0 * reshape(GsPRS(ii,:), [3,1]);
        GpRad(ii,:) = cross(GrRad(ii,:), GsRad(ii,:));
        
        GrRad(ii,:) = GrRad(ii,:) / norm(GrRad(ii,:));
        GsRad(ii,:) = GsRad(ii,:) / norm(GsRad(ii,:));
        GpRad(ii,:) = GpRad(ii,:) / norm(GpRad(ii,:));
        
        [GsRad(ii,:)] = inPlaneRot(GsRad(ii,:), GrRad(ii,:), theta(ii));
        
        R(ii, :, 2) = [GrRad(ii,1); -GrRad(ii,2); -GrRad(ii,3)];
        [dSliceNormalSag, dSliceNormalCor, dSliceNormalTra] = checkNormVec(GsRad(ii,1), GsRad(ii,2), GsRad(ii,3));
        R(ii, :, 3) = [dSliceNormalSag; -dSliceNormalCor;  -dSliceNormalTra];
        R(ii, :, 1) = cross(R(ii, :, 2), R(ii, :, 3));
    end
end

function [dSliceNormalSag, dSliceNormalCor, dSliceNormalTra] = checkNormVec(dSliceNormalSag, dSliceNormalCor, dSliceNormalTra)
    % meaningless, just follow the checkNormalVector() in Siemens, VE11C
    dLargestComponent = dSliceNormalTra;

    if (abs(dSliceNormalCor) > abs(dLargestComponent))
        dLargestComponent = dSliceNormalCor;
    end

    if (abs (dSliceNormalSag) > abs (dLargestComponent))
        dLargestComponent = dSliceNormalSag;
    end

    if (dLargestComponent < 0.0)    %invalid slice orientation
        dSliceNormalSag = -1.0 * dSliceNormalSag;
        dSliceNormalCor = -1.0 * dSliceNormalCor;
        dSliceNormalTra = -1.0 * dSliceNormalTra;
    end
end