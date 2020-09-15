function quat = rot2quat(rotm)
    if ( (size(rotm,1) ~= 3) || (size(rotm,2) ~= 3) )
        error('rot2quat: %s', WBM.wbmErrorMsg.WRONG_MAT_DIM);
    end
    quat    = zeros(4,1);
    epsilon = 1e-12; % min. value to treat a number as zero ...

    %% Compute the corresponding (unit) quaternion from a given rotation matrix R:
    % The transformation uses the computational efficient algorithm of Stanley.
    % To be numerically robust, the code determines the set with the maximum
    % divisor for the calculation.
    % For further details about the Stanley Algorithm, see:
    %   [1] Optimal Spacecraft Rotational Maneuvers, John L. Junkins & James D. Turner, Elsevier, 1986, pp. 28-29, eq. (2.57)-(2.59).
    %   [2] Theory of Applied Robotics: Kinematics, Dynamics, and Control, Reza N. Jazar, 2nd Edition, Springer, 2010, p. 110, eq. (3.149)-(3.152).
    % Note: There exist also an optimized version of the Stanley method and is the fastest
    %       possible computation method for Matlab, but it does not cover all special cases.
    % Further details about the fast calculation can be found at:
    %   [3] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       p. 36, formula (2.30).
    tr = rotm(1,1) + rotm(2,2) + rotm(3,3);
    if (tr > epsilon) % if tr > 0:
            % scalar part:
            quat(1,1) = 0.5*sqrt(tr + 1);
            s_inv = 1/(quat(1,1)*4);
            % vector part:
            quat(2,1) = (rotm(3,2) - rotm(2,3))*s_inv;
            quat(3,1) = (rotm(1,3) - rotm(3,1))*s_inv;
            quat(4,1) = (rotm(2,1) - rotm(1,2))*s_inv;
    else % if tr <= 0, find the greatest diagonal element for calculating
         % the scale factor s and the vector part of the quaternion:
        if ( (rotm(1,1) > rotm(2,2)) && (rotm(1,1) > rotm(3,3)) )
            quat(2,1) = 0.5*sqrt(rotm(1,1) - rotm(2,2) - rotm(3,3) + 1);
            s_inv = 1/(quat(2,1)*4);

            quat(1,1) = (rotm(3,2) + rotm(2,3))*s_inv;
            quat(3,1) = (rotm(2,1) + rotm(1,2))*s_inv;
            quat(4,1) = (rotm(3,1) + rotm(1,3))*s_inv;
        elseif (rotm(2,2) > rotm(3,3))
            quat(3,1) = 0.5*sqrt(rotm(2,2) - rotm(3,3) - rotm(1,1) + 1);
            s_inv = 1/(quat(3,1)*4);

            quat(1,1) = (rotm(1,3) - rotm(3,1))*s_inv;
            quat(2,1) = (rotm(2,1) + rotm(1,2))*s_inv;
            quat(4,1) = (rotm(3,2) + rotm(2,3))*s_inv;
        else
            quat(4,1) = 0.5*sqrt(rotm(3,3) - rotm(1,1) - rotm(2,2) + 1);
            s_inv = 1/(quat(4,1)*4);

            quat(1,1) = (rotm(2,1) - rotm(1,2))*s_inv;
            quat(2,1) = (rotm(3,1) + rotm(1,3))*s_inv;
            quat(3,1) = (rotm(3,2) + rotm(2,3))*s_inv;
        end
    end
end
