function [mesh, beam] = preprocessingBCs(mesh, beam)

% NUMBER OF CONSTRAINED DOFS DUE TO DIRICHLET BCs 
% (for standard case of 1-beam structures)

mesh.fixedDofsInds = cell(2,1);
mesh.numFixedDofs = 0;

mesh.removedLamuInds = cell(2,1);
mesh.numRemovedLamu = 0;

if length(beam.BCs) == 2
    nodeIds = [1; mesh.numNodes];
    for i = 1:2        
        if strcmp(beam.BCs{i}, 'clamped')
            localDofsIds = 1:12;
            localLamuIds = 1:6;
        elseif strcmp(beam.BCs{i}, 'fixed')
            localDofsIds = [1,2,3];
            localLamuIds = [];
        elseif strcmp(beam.BCs{i}, 'free')
            localDofsIds = [];
            localLamuIds=[];
        elseif strcmp(beam.BCs{i}, 'custom')
            if isempty(beam.customFixedDofs{i})
                if i == 1
                    tag = 'left';
                else
                    tag = 'right';
                end
                fprintf(['\nPlease provide an array of local dof number of fixed dofs at the ',tag,' end.\n']);
                fprintf('\nTreated as free-end. The system matrix might be singular.\n');
            else
                localDofsIds = beam.customFixedDofs{i};                
                % fixed all d1, d2, d3
                if all(ismember(4:12,beam.customFixedDofs{i}))
                    localLamuIds=1:6;
                % fixed d1 and d2
                elseif all(ismember(4:9,beam.customFixedDofs{i}))
                    localLamuIds=[1,2,6];
                % fixed d1 and d3
                elseif all(ismember([4,5,6,10,11,12],beam.customFixedDofs{i}))
                    localLamuIds=[1,3,5];
                % fixed d2 and d3
                elseif all(ismember(7:12,beam.customFixedDofs{i}))
                    localLamuIds=[2,3,4];
                % fixed d1 only
                elseif all(ismember([4,5,6],beam.customFixedDofs{i}))
                    localLamuIds = 1;
                % fixed d2 only
                elseif all(ismember([7,8,9],beam.customFixedDofs{i}))
                    localLamuIds = 2;
                % fixed d3 only
                elseif all(ismember([10,11,12],beam.customFixedDofs{i}))
                    localLamuIds = 3;
                end
            end
        else
            if i == 1
                tag = 'left';
            else
                tag = 'right';
            end
            fprintf(['\nThe provided boundary condition at the ',tag,' end is not supported yet.\n']);
            fprintf('\nTreated as free-end. The system matrix might be singular.\n');
            localDofsIds = [];
        end

        % fixed dofs (primal fields)
        mesh.fixedDofsInds{i} = mesh.ID(localDofsIds,nodeIds(i));
        mesh.numFixedDofs = mesh.numFixedDofs + length(mesh.fixedDofsInds{i});
        % irrelevant Lagrange multipliers
        mesh.removedLamuInds{i} = mesh.IDlamu(localLamuIds,nodeIds(i));
        mesh.numRemovedLamu = mesh.numRemovedLamu + length(mesh.removedLamuInds{i});
    end
else
    fprintf('\nPlease provide boundary conditions of the left and right ends in the correct format.\n');
    fprintf('\nNo Dirichlet BCs are enforced (free-free beam). The system matrix might be singular.\n');
end

mesh.numActiveDofs = mesh.ndofsTotal - mesh.numFixedDofs - mesh.numRemovedLamu;