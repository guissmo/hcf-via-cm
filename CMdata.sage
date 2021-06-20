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

def sanityCheck_smoothden(KrPoly, factorbound=10^6, smallfaclimit=None):
    den = KrPoly.denominator();
    F = factorupto(den, bound=factorbound);
    if den == 1:
        return ({1}, 1)
    return (F.keys(), sorted(F.keys())[-1])
    
def sagePolyToPariPoly(poly, sageprelt_to_pari):
  polstr = str(poly).replace("y","x").replace("alphar","("+str(sageprelt_to_pari)+")")
  pol = gp(polstr)
  return(pol)


def pariPolyToSagePoly(poly, pariprelt_to_sage):
  x = PolynomialRing(pariprelt_to_sage.parent(),'x').gen()
  y = pariprelt_to_sage;
  n = Integer(gp.poldegree(poly))
  r = 0
  for i in range(0,n+1):
    coef = gp.polcoef(poly, i)
    m = Integer(gp.poldegree(coef))
    for j in range(0,m+1):
        coef2 = gp.polcoef(coef, j)
        r += y^j*QQ(coef2)*x^i
  return(r)

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
            strcomplexpoly = str(complexpoly);
            if len(strcomplexpoly) > 50:
                strcomplexpoly_a = strcomplexpoly[0:20]
                strcomplexpoly_b = strcomplexpoly[-20:]
            print("Attempting to recognize "+"".join([strcomplexpoly_a, "...", strcomplexpoly_b])+" as a polynomial in "+str(self.Kr)+"\n");
        pol = complexPolytoKrPoly(self.Kr, complexpoly, complexpoly.coefficients()[0].prec());
        if not pol.is_irreducible():
            # if auto:
            #     return("Reducible")
            raise ValueError(f"Polynomial is reducible.")
        return(pol);

    def verify_polynomial(self, pol, factorbound=1000000, smallfaclimit=10000, checkgaloisuntil=100, check_conductor=False, check_conductor_until=1000, auto=False):
        tmp = sanityCheck_smoothden(pol, factorbound=factorbound, smallfaclimit=smallfaclimit)
        pol_gp = sagePolyToPariPoly(pol, self.sageprelt_gp)
        if len(tmp) == 2:
            (F, maxfac) = tmp;
            if maxfac > smallfaclimit:
                if auto:
                    return("Contains large factor in denominator.")
                raise ValueError(f"Large factor found enough! Largest factor was {maxfac} > {smallfaclimit}.")
            probgal = sanityCheck_probsgal(self, pol_gp, F=F, ub=checkgaloisuntil)
            if probgal != 1:
                if auto:
                    return("Is probably not Galois")
                raise ValueError(f"Probably not Galois. :-(")
        if check_conductor:
            try:
                cond = rnfConductor(self, pol_gp, tulong=check_conductor_until)
            except Exception as e:
                raise NotImplementedError(f"Polynomial is not supported by PARI.\nError message: {e}")
            return cond
        return True

    def polredbest(self, pol, convert_to_gp=False):
        pol_gp = sagePolyToPariPoly(pol, self.sageprelt_gp)
        polred_gp = gp.lift(gp.rnfpolredbest(self.Kr_gp, pol_gp))
        ret = pariPolyToSagePoly(polred_gp, self.pariprelt)
        return(ret)

    def find_defining_polynomial_and_verify(self, invs, bp=None, prec=100, verbose=False, factorbound=1000000, smallfaclimit=10000, checkgaloisuntil=100, check_conductor=False, check_conductor_until=1000, autoretry=0):
        r"""find_defining_polynomial_and_verify
        Find the Hilbert class polynomial of the quadratic CM field self.Kr

        + COMPUTING PHASE
          - invs: the invariants to be used, typically an element of rosenhain_invariants(2)
          - bp (basepoint): the chosen period matrix  to be acted upon by Shimura reciprocity
          - prec: the precision in which to evaluate invs

        + VERIFYING PHASE
          - factorbound: up to which integer must we try to factor the denominator of the polynomial we found
          - smallfaclimit (must be <= factorbound): if a factor of the denominator exceeds this, throw error and say that we found a big factor
          - check_conductor: whether or not to run rnfconductor to verify conductor (slow)

          - autoretry: how many times I have to retry by doubling the precision before giving up
          - verbose: whether or not I should print progress
        """
        pol = self.defining_polynomial(invs, bp=bp, prec=prec, verbose=verbose)
        ver = self.verify_polynomial(pol, factorbound=factorbound, smallfaclimit=smallfaclimit, checkgaloisuntil=checkgaloisuntil, check_conductor=check_conductor, check_conductor_until=check_conductor_until, auto=True)
        if ver in ["Contains large factor in denominator.", "Is probably not Galois"]:
            prec2 = prec*2;
            newsmallfaclimit = 1000000;
            if autoretry:
                if verbose:
                    print(f"! Polynomial obtained {ver.lower()}.")
                    print(f"! Increasing precision from {prec} to {prec2}.")
                    print(f"! Increasing smallfaclimit to 1000000.")
                return(self.find_defining_polynomial_and_verify(invs, bp=None, prec=prec2, verbose=verbose, factorbound=factorbound, smallfaclimit=newsmallfaclimit, checkgaloisuntil=checkgaloisuntil, check_conductor=check_conductor, check_conductor_until=check_conductor_until, autoretry=autoretry-1))
            else:
                raise ValueError("Given up.")
        return(pol, ver)

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