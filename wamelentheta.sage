load("/mnt/data/backups/vc/master/phdthesis/code/wamtheta/wamtheta.sage")

CC = magma("ComplexField(100)");
I = magma.Name(CC, 1);
CC.AssignNames(["I"]);

def sageCCvec_from_magmaCCvec(vec, prec=2000, magmavec=False):
    if magmavec:
        n = Integer(magma.Ncols(vec));
        return( [ComplexField(prec)(vec[i]) for i in range(1,n+1)] );
    return( [ ComplexField(prec)(i) for i in vec ] );

def magmaCC_withprec(prec):
    return( magma("ComplexField("+str(prec)+")") );

def t1t2t3_from_PM(PM, prec=None):
    ret = [PM[0,0], PM[1,1], PM[0,1]];
    if prec == None:
        return( ret );
    return( sageCCvec_from_magmaCCvec(ret, prec=prec) );

def magmaPM_from_sagePM(PM, prec=2000):
    [t1, t2, t3] = t1t2t3_from_PM(PM, prec=prec);
    return( magma.Matrix(magmaCC_withprec(prec), 2, 2, [t1,t3,t3,t2]) );

def magmaZvec_to_magmaCvec(a1a2b1b2, tau, prec=2000):
    
    [z1, z2] = magma.function_call('Zvec_to_Cvec', [a1a2b1b2, tau], nvals=2);
    z = magma.Matrix(magmaCC_withprec(prec), 2, 1, [z1, z2]);
    return(z);

# tau = magmaPM_from_sagePM(Zs[bpind]);

#threetors_zvec = magma.Vector([0,0,1/3,0]);
#z = magmaZvec_to_magmaCvec(threetors_zvec, tau);

def funsqconsq_from_zPM(z, sagePM, prec=2000):
    prechi = prec*2; # maybe should be a bit higher?
    magmaPM = magmaPM_from_sagePM(sagePM, prec=prechi);
    if len(z) == 4:
        z = magmaZvec_to_magmaCvec(z, magmaPM, prec=prechi);
    tau = magmaPM_from_sagePM(sagePM);
    [fs, cs] = magma.function_call('CalculThetas', [z, tau/2], nvals=2);
    fs = [f0,f1,f2,f3] = sageCCvec_from_magmaCCvec(magma.Vector(fs), prec=prec, magmavec=True);
    cs = [c0,c1,c2,c3] = sageCCvec_from_magmaCCvec(magma.Vector(cs), prec=prec, magmavec=True);
    fun = [
        (f0*c0 + f1*c1 + f2*c2 + f3*c3)/4,
        (f0*c1 + f1*c0 + f2*c3 + f3*c2)/4,
        (f0*c2 + f1*c3 + f2*c0 + f3*c1)/4,
        (f0*c3 + f1*c2 + f2*c1 + f3*c0)/4,
        (f0*c0 - f1*c1 + f2*c2 - f3*c3)/4,
        (f0*c1 - f1*c0 + f2*c3 - f3*c2)/4,
        (f0*c2 - f1*c3 + f2*c0 - f3*c1)/4,
        (f0*c3 - f1*c2 + f2*c1 - f3*c0)/4,
        (f0*c0 + f1*c1 - f2*c2 - f3*c3)/4,
        (f0*c1 + f1*c0 - f2*c3 - f3*c2)/4,
        (f0*c2 + f1*c3 - f2*c0 - f3*c1)/4,
        (f0*c3 + f1*c2 - f2*c1 - f3*c0)/4,
        (f0*c0 - f1*c1 - f2*c2 + f3*c3)/4,
        (f0*c1 - f1*c0 - f2*c3 + f3*c2)/4,
        (f0*c2 - f1*c3 - f2*c0 + f3*c1)/4,
        (f0*c3 - f1*c2 - f2*c1 + f3*c0)/4
    ];
    con = [
        (c0*c0 + c1*c1 + c2*c2 + c3*c3)/4,
        (c0*c1 + c1*c0 + c2*c3 + c3*c2)/4,
        (c0*c2 + c1*c3 + c2*c0 + c3*c1)/4,
        (c0*c3 + c1*c2 + c2*c1 + c3*c0)/4,
        (c0*c0 - c1*c1 + c2*c2 - c3*c3)/4,
        (c0*c1 - c1*c0 + c2*c3 - c3*c2)/4,
        (c0*c2 - c1*c3 + c2*c0 - c3*c1)/4,
        (c0*c3 - c1*c2 + c2*c1 - c3*c0)/4,
        (c0*c0 + c1*c1 - c2*c2 - c3*c3)/4,
        (c0*c1 + c1*c0 - c2*c3 - c3*c2)/4,
        (c0*c2 + c1*c3 - c2*c0 - c3*c1)/4,
        (c0*c3 + c1*c2 - c2*c1 - c3*c0)/4,
        (c0*c0 - c1*c1 - c2*c2 + c3*c3)/4,
        (c0*c1 - c1*c0 - c2*c3 + c3*c2)/4,
        (c0*c2 - c1*c3 - c2*c0 + c3*c1)/4,
        (c0*c3 - c1*c2 - c2*c1 + c3*c0)/4
    ];
    return([fun, con, fs, cs])

def upolvpolsq_from_zPM(z, sagePM, prec=2000):

    [fun, con, fs, cs] = funsqconsq_from_zPM(z, sagePM, prec=prec);

    [t0_z_L_sq, t1_z_L_sq, t2_z_L_sq, t3_z_L_sq, t4_z_L_sq, t5_z_L_sq, t6_z_L_sq, t7_z_L_sq, t8_z_L_sq, t9_z_L_sq, t10_z_L_sq, t11_z_L_sq, t12_z_L_sq, t13_z_L_sq, t14_z_L_sq, t15_z_L_sq] = fun;
    [t0_0_L_sq, t1_0_L_sq, t2_0_L_sq, t3_0_L_sq, t4_0_L_sq, t5_0_L_sq, t6_0_L_sq, t7_0_L_sq, t8_0_L_sq, t9_0_L_sq, t10_0_L_sq, t11_0_L_sq, t12_0_L_sq, t13_0_L_sq, t14_0_L_sq, t15_0_L_sq] = con;

    u1 = (t0_0_L_sq * t2_0_L_sq * t8_0_L_sq) / (t4_0_L_sq * t6_0_L_sq * t12_0_L_sq) * ( t10_z_L_sq^2 / t14_z_L_sq^2 )
    u2 = (t1_0_L_sq * t3_0_L_sq * t9_0_L_sq) / (t4_0_L_sq * t6_0_L_sq * t12_0_L_sq) * ( t11_z_L_sq^2 / t14_z_L_sq^2 )
    upol = x*(x-1) - (x - 1)*u1 + x*u2

    [x2, w2, y2, z2] = fs;
    [a2, b2, c2, d2] = cs;

    Ap = a2 + b2 + c2 + d2
    Bp = a2 + b2 - c2 - d2
    Cp = a2 - b2 + c2 - d2
    Dp = a2 - b2 - c2 + d2
    Ep = Ap*Bp*Cp*Dp / ( (a2*d2 - b2*c2)*(a2*c2 - b2*d2)*(a2*b2 - c2*d2) )
    F = (a2^2 - b2^2 - c2^2 + d2^2) / ( a2*d2 - b2*c2 )
    G = (a2^2 - b2^2 + c2^2 - d2^2) / ( a2*c2 - b2*d2 )
    H = (a2^2 + b2^2 - c2^2 - d2^2) / ( a2*b2 - c2*d2 )
    xyztabcd = ( F*(x2*w2 + y2*z2) + G*(x2*z2 + y2*w2) + H*(x2*y2 + z2*w2) - (x2^2 + y2^2 + z2^2 + w2^2) ) / (2 * Ep)

    tmp37111504812 = t3_z_L_sq * t15_z_L_sq * t0_0_L_sq * t12_0_L_sq - t0_z_L_sq * t1_z_L_sq * t0_0_L_sq * t1_0_L_sq + xyztabcd;
    tmp27101514912 = t2_z_L_sq * t15_z_L_sq * t1_0_L_sq * t12_0_L_sq - t1_z_L_sq * t2_z_L_sq * t1_0_L_sq * t2_0_L_sq + xyztabcd;
    tmp2310110189 = t2_z_L_sq * t3_z_L_sq * t0_0_L_sq * t1_0_L_sq - xyztabcd;

    Y12Y13 = (t0_0_L_sq^2 * t1_0_L_sq^2 * t2_0_L_sq^2 * t3_0_L_sq * t8_0_L_sq^2 * t9_0_L_sq^2 * t15_0_L_sq) / (t4_0_L_sq^5 * t6_0_L_sq^4 * t12_0_L_sq^5) * (t10_z_L_sq / t14_z_L_sq^3) * tmp37111504812;
    Y12Y23 = (t0_0_L_sq^2 * t1_0_L_sq^2 * t2_0_L_sq * t3_0_L_sq^2 * t8_0_L_sq^2 * t9_0_L_sq^2 * t15_0_L_sq) / (t4_0_L_sq^5 * t6_0_L_sq^4 * t12_0_L_sq^5) * (t11_z_L_sq / t14_z_L_sq^3) * tmp27101514912;
    Y13Y23 = (t0_0_L_sq^2 * t1_0_L_sq^2 * t2_0_L_sq * t3_0_L_sq * t8_0_L_sq^2 * t9_0_L_sq^2 * t15_0_L_sq^2) / (t4_0_L_sq^5 * t6_0_L_sq^4 * t12_0_L_sq^5) * (t7_z_L_sq / t14_z_L_sq^3) * tmp2310110189;

    Y12sq = (t0_0_L_sq^2 * t1_0_L_sq^2 * t2_0_L_sq^2 * t3_0_L_sq^2 * t8_0_L_sq^2 * t9_0_L_sq^2) / (t4_0_L_sq^4 * t6_0_L_sq^4 * t12_0_L_sq^4)  *  (t10_z_L_sq * t11_z_L_sq * t15_z_L_sq) / (t14_z_L_sq^3);
    Y13sq = (t0_0_L_sq^3 * t1_0_L_sq^2 * t2_0_L_sq^2 * t8_0_L_sq^3 * t9_0_L_sq^2 * t15_0_L_sq^2) / (t4_0_L_sq^5 * t6_0_L_sq^4 * t12_0_L_sq^5) * (t3_z_L_sq * t7_z_L_sq * t10_z_L_sq) / (t14_z_L_sq^3);
    Y23sq = (t0_0_L_sq^2 * t1_0_L_sq^3 * t3_0_L_sq^2 * t8_0_L_sq^2 * t9_0_L_sq^3 * t15_0_L_sq^2) / (t4_0_L_sq^5 * t6_0_L_sq^4 * t12_0_L_sq^5)  *  (t2_z_L_sq * t7_z_L_sq * t11_z_L_sq) / (t14_z_L_sq^3);

    lam = (t0_0_L_sq * t1_0_L_sq) / (t2_0_L_sq * t3_0_L_sq);

    v1sq = 1/(lam - 1)^2 * (Y12sq + Y13sq + 2*Y12sq*Y13sq);
    v2sq = 1/(lam - 0)^2 * (Y12sq + Y23sq + 2*Y12sq*Y23sq);
    v1v2 = 1/(lam*(lam-1)) * (Y12sq - Y12Y13 - Y12Y13 + Y13Y23);

    vpolsq  =  (x - 1)^2 * v1sq  +  x * v2sq  +  x*(x-1) * v1v2;

    return([upol, vpolsq]);

[upol, vpolsq] = upolvpolsq_from_zPM([0,0,1/3,0], Zs[bpind], prec=1000)

def u0_from_upol(upol):
    return( upol(x=0) );

def u0_from_zPM(z, sagePM, prec=2000):
    [upol, vpolsq] = upolvpolsq_from_zPM(z, sagePM, prec=prec);
    return( u0_from_upol(upol) );

u0_from_zPM([0,0,1/3,0], Zs[bpind], prec=92) 
