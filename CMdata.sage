gp.read(FGAGDIR+"/fgag.gp");
gp.read(FGAGDIR+"/fgagshimuray.gp");

def ShGrpKr_from_Vm(V, m=None):
    ShGrp = gp.fgaginitshimuray(V, m, 1);
    Containment = gp.cmkcontainshcf(ShGrp, m, 1);
    if Containment == 0:
        raise ValueError("HK(1) does not contain CMK(2)");
    Kr = (gp.shraycontextcm(ShGrp))[7];
    gp.set('ShGrp', str(ShGrp));
    gp.set('Kr', str(Kr));
    return ([ShGrp, Kr]);

def sagepariPrimitiveElement(pariKr, pariprelt, sageKr, sageprelt):
  transfor = gp.Vec( gp.nfisisom( pariKr, gp.bnfinit( sageKr.polynomial() ) )[2] )
  fortrans = gp.Vec( gp.nfisisom( gp.bnfinit( sageKr.polynomial() ), pariKr )[2] )
  pariprelt_to_sage = sum([ Rational(transfor[len(transfor)-i])*sageprelt^i for i in range(0, len(transfor)) ]);
  sageprelt_to_pari = sum([ Rational(fortrans[len(fortrans)-i])*pariprelt^i for i in range(0, len(fortrans)) ]);
  return [pariprelt_to_sage, sageprelt_to_pari]

def convertIdeal(pariKr, sageKr, pariprelt_to_sage, idl):
  idl = gp.idealmoddivisor( gp.bnrinit(pariKr, 2), idl );
  [ a, fraka ] = gp.idealtwoelt( pariKr, idl )
  fraka_powerbasis = gp.Vecrev( (pariKr[7][7]/pariKr[7][7][1]) * fraka )
  ret = sageKr.ideal ( Rational(a), sum([ Rational(fraka_powerbasis[i+1])*pariprelt_to_sage^i for i in range(0, len(fraka_powerbasis)) ]) )
  return(ret);

# amIreal, used by rootstoPoly's reality check
def amIreal(pol):
    imag = pol.map_coefficients(imag_part)
    lst = list(imag)
    for coef in imag:
        if abs(coef) > 10^-10:
            return False
        else:
            return True

# rootstoPoly, appears in 3.1
def rootstoPoly(rts, realitycheck=False):
    x = polygen(rts[0].parent())
    prod = 1
    for i in range(0,len(rts)):
        prod = prod * (x - rts[i])
    if realitycheck:
        if amIreal(prod):
             return prod.map_coefficients(real_part)
        else:
             print(prod);
             raise ValueError("The polynomial does not seem to have real coefficients, but it should!")
    return(prod);

def complexPolytoKrPoly(Kr, pol, pr):
    CCpr = ComplexField(pr)
    emb = Kr.embedding(CCpr)
    newpol = recognize_polynomial(pol, Kr, emb=emb)
    return newpol

def factorupto(N, bound=10^6):
    F_gp = gp.factor(N, bound);
    [row, col] = gp.matsize(F_gp);
    F = {};
    for i in range(0,Integer(row)):
        F[F_gp[i+1,1]] = F_gp[i+1,2]
    return F;

def sanityCheck_smoothden(KrPoly, factorbound=10^6, smoothbound=None):
    den = KrPoly.denominator();
    F = factorupto(den, bound=factorbound);
    if den == 1:
        return True
    return (F.keys(), sorted(F.keys())[-1])

def sagePolyToPariPoly(poly, sageprelt_to_pari):
  polstr = str(poly).replace("y","x").replace("alphar","("+str(sageprelt_to_pari)+")")
  pol = gp(polstr)
  return(pol)

def sanityCheck_probsgal(sgdata, poly, F=None, ub=100):
    if F is None:
        (F, _) = sanityCheck_smoothden(poly);
    if not isinstance(poly, sage.interfaces.gp.GpElement): #then it is a SAGE poly
        poly = sagePolyToPariPoly(poly, sgdata.sageprelt_gp)
    gp_F = gp( [ i for i in F] )
    return( gp.isprobablygalois(sgdata.Kr_gp, poly, gp_F, ub) );

def rnfConductor(sgdata, poly, tulong=10000):
  cond = gp.rnfconductor(sgdata.Kr_gp, [poly, tulong])
  if cond == 0:
    raise ValueError;
    return(0)
  return(cond[1])

####

class ShGpData():

    def __init__(self, V, m=2, bpind=None, alphar=true):

        self.V = V;
        self.m = m;

        [self.ShGrp_gp, self.Kr_gp] = ShGrpKr_from_Vm(self.V,m=self.m);
        
        self.K = CM_Field(V);
        self.Zs = Zs = list(self.K.period_matrices_iter());

        while true:
            self.randombp = false;
            if bpind == None:
                self.randombp = true;
                self.bpind = randrange(0, len(self.Zs))
            else:
                self.bpind = bpind;
            self.bp = self.Zs[self.bpind];
            self.cmt = self.bp.CM_type();
            self.Kr = self.cmt.reflex_field();
            if not alphar:
                break
            else:
                if self.Kr.variable_name() == 'alphar':
                    break

        self.yr = self.Kr.absolute_generator();
        self.yr_gp = gp("y");

        [self.pariprelt, self.sageprelt_gp] = sagepariPrimitiveElement(self.Kr_gp, self.yr_gp, self.Kr, self.yr);

        self.cosets_gp = gp.shraycosetsmacronaive(self.ShGrp_gp, self.m);

        self.cosets = [];

        for coset_gp in self.cosets_gp:
            self.coset = [];
            coset_dlog_gp = coset_gp[1];
            for dlgidl in coset_gp[2]:
                self.coset.append( ShGpCosetRep( \
                   convertIdeal(self.Kr_gp, self.Kr, self.pariprelt, dlgidl[2]), \
                   [ZZ(gp.lift(z)) for z in dlgidl[1]]
                ));
            self.cosets.append( ShGpCoset(self.coset, coset_dlog_gp) );

    def complex_defining_polynomial(self, invs, bp=None, prec=100, verbose=False):
        if bp == None:
            bp = self.Zs[self.bpind];
        if verbose:
            print("We are trying to compute a degree "+str(len(self.cosets))+" polynomial.\n");
        roots = [coset.sum_of_invs(invs, bp, prec=prec, verbose=verbose) for coset in self.cosets];
        return ( rootstoPoly(roots) );

    def defining_polynomial(self, invs, bp=None, prec=100, verbose=False):
        complexpoly = self.complex_defining_polynomial(invs, bp, prec=prec, verbose=verbose);
        if verbose:
            print("Attempting to recognize "+str(complexpoly)+" as a polynomial in "+str(self.Kr)+"\n");
        return(complexPolytoKrPoly(self.Kr, complexpoly, complexpoly.coefficients()[0].prec()));

    def verify_polynomial(self, pol, factorbound=1000000, smoothbound=10000, checkgaloisuntil=100, check_conductor=False, check_conductor_until=1000):
        tmp = sanityCheck_smoothden(pol, factorbound=factorbound, smoothbound=smoothbound);
        if len(tmp) != 2:
            return tmp
        (F, maxfac) = tmp;
        if maxfac > smoothbound:
            raise ValueError(f"Not smooth enough! Largest factor was {maxfac} > {smoothbound}.")
        pol_gp = sagePolyToPariPoly(pol, self.sageprelt_gp)
        probgal = sanityCheck_probsgal(self, pol_gp, F=F, ub=checkgaloisuntil)
        if probgal != 1:
            raise ValueError(f"Probably not Galois. :-(")
        if check_conductor:
            cond = rnfConductor(self, pol_gp, tulong=check_conductor_until)
        print(cond)
        return True

class ShGpCoset():

    def __init__(self, reps, coset_dlog_gp):
        self.reps = reps;
        self.coset_dlog = [ZZ(gp.lift(z)) for z in coset_dlog_gp];

    def __repr__(self):
        return( str("Coset with dlog "+str(self.coset_dlog))+" containing "+str(len(self.reps))+" elements" );

    def sum_of_invs(self, inv, bp, n=8, prec=100, prec2=None, verbose=verbose):
        ret = 0;
        if verbose:
            print("  Going through elements of "+self.__repr__()+"\n");
        for rep in self.reps:
            [rep.U, rep.m, rep.u] = bp.Shimura_reciprocity(rep.idl, n=n, period_matrix=True);
        for rep in self.reps:
            ret = ret+rep.evaluate_inv(inv, prec=prec, verbose=verbose);
        return(ret);

class ShGpCosetRep():

    def __init__(self, idl, repdlog):
        self.idl = idl;
        self.repdlog = repdlog;
        self.U = None;
        self.m = None;
        self.u = None;
        self.theta = None;
    
    def __repr__(self):
        return( "(" + str(self.idl.gens()[0]) + ", ...) " + str(self.repdlog) );

    def evaluate_inv(self, inv, prec=None, verbose=verbose):
        if verbose:
            print("    Evaluating "+str(inv^self.u)+" with precision "+str(prec)+" with z = 0 and tau = "+str(self.U)+"\n");
        self.theta = (inv^self.u)(self.U, prec=prec, use_magma=True);
        return(self.theta);