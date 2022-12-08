function rotated_image = imrotatepad(image, rot_angle_degree)
    RA = imref2d(size(image));    
    tform = affine2d([cosd(rot_angle_degree)    -sind(rot_angle_degree)     0; ...
                      sind(rot_angle_degree)     cosd(rot_angle_degree)     0; ...
                      0                          0                          1]);
      Rout = images.spatialref.internal.applyGeometricTransformToSpatialRef(RA,tform);
      Rout.ImageSize = RA.ImageSize;
      xTrans = mean(Rout.XWorldLimits) - mean(RA.XWorldLimits);
      yTrans = mean(Rout.YWorldLimits) - mean(RA.YWorldLimits);
      Rout.XWorldLimits = RA.XWorldLimits+xTrans;
      Rout.YWorldLimits = RA.YWorldLimits+yTrans;
      rotated_image = imwarp(image, tform, 'OutputView', Rout, 'interp', 'cubic', 'fillvalues', 255);
  end