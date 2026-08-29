function answer=isworkernode()
% BasisREMY Octave shim shadowing Spinach's kernel/utilities/isworkernode.m,
% which asks MATLAB's Parallel Computing Toolbox. Answering "worker node"
% makes Spinach skip every pool / parallel branch (create, hamiltonian,
% evolution, relaxation, basis) and run serially, which is what Octave does.
answer=true;
end
