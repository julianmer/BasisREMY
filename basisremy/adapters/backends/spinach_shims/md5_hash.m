function hashstr=md5_hash(A)
% BasisREMY Octave shim shadowing Spinach's kernel/utilities/md5_hash.m, which
% is built on MATLAB p-code (serializeToBytes / digestMD5) that Octave cannot
% load. Spinach only needs a stable digest of numeric / char / cell / struct
% inputs for its cache keys, so a plain serialisation hashed with hash() does.
hashstr=hash('md5',char(ser(A)));
end
function b=ser(A)
if isnumeric(A)||islogical(A)
    b=[uint8(class(A)) uint8(0) typecast(double(size(A)),'uint8') ...
       typecast(double(real(full(A(:)))).','uint8') typecast(double(imag(full(A(:)))).','uint8')];
elseif ischar(A)
    b=[uint8('char') uint8(0) typecast(double(size(A)),'uint8') uint8(A(:).')];
elseif iscell(A)
    b=[uint8('cell') uint8(0) typecast(double(size(A)),'uint8')];
    for k=1:numel(A), b=[b ser(A{k})]; end
elseif isstruct(A)
    f=sort(fieldnames(A)); b=[uint8('struct') uint8(0) typecast(double(size(A)),'uint8')];
    for k=1:numel(A), for j=1:numel(f), b=[b uint8(f{j}) uint8(0) ser(A(k).(f{j}))]; end; end
elseif isa(A,'function_handle')
    b=[uint8('fh') uint8(0) uint8(func2str(A))];
else
    b=[uint8(class(A)) uint8(0) uint8(disp(A))];
end
end
