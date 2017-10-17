% Dielectric tensor definition function
% Last edit: 10-12-2017
% Editor: Casey Icenhour (NCSU)

function tensor = dielectricTensor(type,const)

    if strcmp(type,'vacuum')
        tensor = eye(3);
    elseif strcmp(type,'constant')
        tensor = const*eye(3);
    else
        errmsg('The only currently valid tensor types are vacuum and constant!')
    end

end

