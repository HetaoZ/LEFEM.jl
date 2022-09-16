function elast_matrix_1d(E, C, ptype)
    if ptype in ("bar" , "beam")
        D = E * C
    else
        error("undef")
    end
    return D
end

elast_matrix_bar(E, A) = E*A

elast_matrix_beam(E, I) = E*I

function elast_matrix_2d(E, ν, ptype)
    if ptype == "pstrain"
        D = elast_matrix_pstrain(E, ν)
    elseif ptype == "pstress"
        D = elast_matrix_pstress(E, ν)
    elseif ptype == "sym"
        D = elast_matrix_sym(E, ν)
    else
        error("undef")
    end
    return D
end

function elast_matrix_pstrain(E, ν)
    E0 = E / (1 - ν^2)
    ν0 = ν/ (1 - ν)
    return elast_matrix_pstress(E0, ν0)
end

function elast_matrix_pstress(E, ν)
    D0 = E / (1 - ν^2)
    D = D0 * [1 ν 0; ν 1 0; 0 0 (1-ν)/2]
    return D
end

function elast_matrix_sym(E, ν)
    D0 = E*(1-ν) / ((1+ν)*(1-2*ν))
    a = ν / (1-ν)
    b = (1-2*ν)/(2*(1-ν))
    D = D0 * 
    [1 a 0 a; 
    a 1 0 a; 
    0 0 b 0; 
    a a 0 1]
    return D
end

function elast_matrix_3d(E, ν)
    D0 = E*(1-ν) / ((1+ν)*(1-2*ν))
    a = ν / (1-ν)
    b = (1-2*ν)/(2*(1-ν))
    D = D0 * 
    [1 a a 0 0 0; 
    a 1 a 0 0 0; 
    a a 1 0 0 0; 
    0 0 0 b 0 0;
    0 0 0 0 b 0;
    0 0 0 0 0 b]
    return D
end
