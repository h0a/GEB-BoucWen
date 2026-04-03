function [Ac, Hc, hc, rqc] = assembleLinTangentJointCons(meshes, beams, q, lambda)

% allocation
Ac = zeros(meshes.ndofsQtot);
Hc = zeros(meshes.numLambda4coupling,meshes.ndofsQtot);

rqc = zeros(meshes.ndofsQtot,1);
hc = zeros(meshes.numLambda4coupling,1);

% find all beams which have slave node(s)
[~,colId] = find(beams.IENcoup < 0);
slaveBeamNo = sort(unique(colId));

clambda2 = 0;
for cb = 1:length(slaveBeamNo)      % loop over all beams which have slave node
    b = slaveBeamNo(cb);
    for j = 1:2                     % loop over 2 ends of the current beam
        if beams.IENcoup(j,b) < 0
            meshSlave = meshes.meshi{b};
            beamSlave = beams.beami{b};

            bmaster_id = beams.IENmaster4coup{j,b}(1);
            master_flag = beams.IENmaster4coup{j,b}(2);

            beamMaster = beams.beami{bmaster_id};
            meshMaster = meshes.meshi{bmaster_id};

            % extract lambda for the current joint
            clambda1 = clambda2 + 1;
            if beams.IENcoup(j,b) == -1                
                clambda2 = clambda2 + 12;
            elseif beams.IENcoup(j,b) == -2
                clambda2 = clambda2 + 3;
            end
            lambda_ms = lambda(clambda1:clambda2);

            % extract q of the current slave node
            cqS = 0;
            if b > 1
                for bb = 1:b-1
                    cqS = cqS + meshes.meshi{bb}.ndofs;
                end
            end
            
            if j == 1
                id_qS = cqS+1:cqS+12;             % 1st node
                id_Qs = meshSlave.ID(:,1);
            else
                cqS = cqS + meshSlave.ndofs;
                id_qS = cqS-11:cqS;               % last node
                id_Qs = meshSlave.ID(:,end);
            end
            qSlave  = q(id_qS);

            % extract q of the current master node
            cqM = 0;
            if bmaster_id > 1
                for bb = 1:bmaster_id-1
                    cqM = cqM + meshes.meshi{bb}.ndofs;
                end
            end

            if master_flag == 0               % 1st node of master beam
                id_qM = cqM+1:cqM+12;
                id_Qm = meshMaster.ID(:,1);
            else                              % last node of master beam
                cqM = cqM + meshMaster.ndofs;
                id_qM = cqM-11:cqM;
                id_Qm = meshMaster.ID(:,end);
            end
            qMaster = q(id_qM);

            % linearization of the coupling constraint
            if beams.IENcoup(j,b) == -1
                Qslave = beamSlave.q(id_Qs,1);
                Qmaster = beamMaster.q(id_Qm,1);

                [Kms, Hmaster, Hslave, hms] = linearizeRigidJointeCons(qMaster, qSlave, Qmaster, Qslave, lambda_ms);
            elseif beams.IENcoup(j,b) == -2
                [Kms, Hmaster, Hslave, hms] = linearizeHingedJointeCons(qMaster, qSlave);
            end
            rm = Hmaster' * lambda_ms;
            rs = Hslave'  * lambda_ms;

            % assign to the global matrix and rhs
            Ac(id_qM,id_qS) = Ac(id_qM,id_qS) + Kms;
            Ac(id_qS,id_qM) = Ac(id_qS,id_qM) + Kms';

            Hc(clambda1:clambda2,id_qS) = Hc(clambda1:clambda2,id_qS) + Hslave;
            Hc(clambda1:clambda2,id_qM) = Hc(clambda1:clambda2,id_qM) + Hmaster;
            
            rqc(id_qS) = rqc(id_qS) + rs;
            rqc(id_qM) = rqc(id_qM) + rm;
            
            hc(clambda1:clambda2) = hc(clambda1:clambda2) + hms;
        end
    end                     % loop over 2 ends of the current beam
end                         % loop over all beams which have slave node
end


function [Kms, Hmaster, Hslave, hms] = linearizeRigidJointeCons(qMaster, qSlave, Qmaster, Qslave, lambda)

    % extract relevant solutions    
    phi01 = qMaster(1:3);
    d11   = qMaster(4:6);
    d21   = qMaster(7:9);
    d31   = qMaster(10:12);
    
    D11   = Qmaster(4:6);
    D21   = Qmaster(7:9);
    D31   = Qmaster(10:12);
    
    phi02 = qSlave(1:3);
    d12   = qSlave(4:6);
    d22   = qSlave(7:9);
    d32   = qSlave(10:12);

    D12   = Qslave(4:6);
    D22   = Qslave(7:9);
    D32   = Qslave(10:12);

    % tangent stiffness matrices
    Kms = [sparse(3,12);
            sparse(3,3) lambda(4) .*speye(3)    lambda(5) .*speye(3)    lambda(6) .*speye(3);
            sparse(3,3) lambda(7) .*speye(3)    lambda(8) .*speye(3)    lambda(9) .*speye(3);
            sparse(3,3) lambda(10).*speye(3)    lambda(11).*speye(3)    lambda(12).*speye(3)];

    % constraint matrices
    Hmaster = [ speye(3)  sparse(3,9);
           sparse(1,3)  d12'    sparse(1,6);
           sparse(1,3)  d22'    sparse(1,6);
           sparse(1,3)  d32'    sparse(1,6);
           sparse(1,6)  d12'    sparse(1,3);
           sparse(1,6)  d22'    sparse(1,3);
           sparse(1,6)  d32'    sparse(1,3);
           sparse(1,9)  d12';
           sparse(1,9)  d22';
           sparse(1,9)  d32'];

    Hslave = [-speye(3)  sparse(3,9);
           sparse(1,3)  d11'    sparse(1,6);
           sparse(1,6)  d11'    sparse(1,3);
           sparse(1,9)  d11';
           sparse(1,3)  d21'    sparse(1,6);
           sparse(1,6)  d21'    sparse(1,3);
           sparse(1,9)  d21';
           sparse(1,3)  d31'    sparse(1,6);
           sparse(1,6)  d31'    sparse(1,3);
           sparse(1,9)  d31'];
    
    % constraint vector
    hms = [phi01-phi02;
          dot(d11,d12) - dot(D11,D12);
          dot(d11,d22) - dot(D11,D22);
          dot(d11,d32) - dot(D11,D32);
          dot(d21,d12) - dot(D21,D12);
          dot(d21,d22) - dot(D21,D22);
          dot(d21,d32) - dot(D21,D32);
          dot(d31,d12) - dot(D31,D12);
          dot(d31,d22) - dot(D31,D22);
          dot(d31,d32) - dot(D31,D32)];
end


function [Kms, Hmaster, Hslave, hms] = linearizeHingedJointeCons(qMaster, qSlave)

    % extract relevant solutions    
    phi01 = qMaster(1:3);    
    phi02 = qSlave(1:3);

    Kms = sparse(12,12);

    % constraint matrices
    Hmaster = [ speye(3)  sparse(3,9)];
    Hslave  = [-speye(3)  sparse(3,9)];

    % constraint vector
    hms = phi01-phi02;
end
