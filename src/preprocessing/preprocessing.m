function [mesh, beam] = preprocessing(mesh, beam)

[mesh, beam] = preprocessingGeneral(mesh, beam);
[mesh, beam] = preprocessingBCs(mesh, beam);