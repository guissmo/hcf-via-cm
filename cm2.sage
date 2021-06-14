# parisetup, appears in 0.1
# - uses the gp function fgaginitshimuray from fgagshimuray
# - uses the gp function cmkcontainshcf from fgagshimuray
# - uses the gp function shraycontextcm from fgagshimuray
def parisetup(V, m):
  ShGrp = gp.fgaginitshimuray(V, m, 1);
  Containment = gp.cmkcontainshcf(ShGrp, m, 1);
  if Containment == 0:
    raise ValueError("HK(1) does not contain CMK(2)");
  Kr = (gp.shraycontextcm(ShGrp))[7];
  gp.set('ShGrp', str(ShGrp));
  gp.set('Kr', str(Kr));
  return ([ShGrp, Kr]);

# sagesetup, appears in 0.2
# see cm2template for output
# - uses several functions from recip
def sagesetup(V, ind=None, Zs=None):
  K = CM_Field(V)
  if Zs == None:
    Zs = list(K.period_matrices_iter())
  bpind = ind;
  if bpind == None:
    bpind = randrange(0, len(Zs))
  bp = Zs[ bpind ]
  CMtype = bp.CM_type()
  Kr = bp.CM_type().reflex_field()
  if Kr.variable_name() == 'alpharp':
    return( sagesetup(V, None, Zs=Zs) )
  print("using bpind %d..." % bpind)
  yr = Kr.absolute_generator()
  return ([K, Zs, [bp, bpind], [Kr, yr]]) ;

# sagepariPrimitiveElement, appears in 0.3
def sagepariPrimitiveElement(pariKr, pariprelt, sageKr, sageprelt):
  transfor = gp.Vec( gp.nfisisom( pariKr, gp.bnfinit( sageKr.polynomial() ) )[2] )
  fortrans = gp.Vec( gp.nfisisom( gp.bnfinit( sageKr.polynomial() ), pariKr )[2] )
  pariprelt_to_sage = sum([ Rational(transfor[len(transfor)-i])*sageprelt^i for i in range(0, len(transfor)) ]);
  sageprelt_to_pari = sum([ Rational(fortrans[len(fortrans)-i])*pariprelt^i for i in range(0, len(fortrans)) ]);
  return [pariprelt_to_sage, sageprelt_to_pari]

# pariCosets, appears in 1.1
def pariCosets(ShGrp, m):
  return( gp.shraycosetsmacronaive(ShGrp, m) );

# cosetLabel, used by sagepariIdealConversion
def cosetLabel(gpvec):
  return( ' '.join([str(gp.lift(g)) for g in gpvec]) );

# convertIdeal, used by sagepariIdealConversion
def convertIdeal(pariKr, sageKr, pariprelt_to_sage, idl):
  idl = gp.idealmoddivisor( gp.bnrinit(pariKr, 2), idl );
  [ a, fraka ] = gp.idealtwoelt( pariKr, idl )
  fraka_powerbasis = gp.Vecrev( (pariKr[7][7]/pariKr[7][7][1]) * fraka )
  ret = sageKr.ideal ( Rational(a), sum([ Rational(fraka_powerbasis[i+1])*pariprelt_to_sage^i for i in range(0, len(fraka_powerbasis)) ]) )
  return(ret);

# sagepariIdealConversion, appears in 1.2
def sagepariIdealConversion(pariKr, sageKr, pariprelt_to_sage, listahan):
  dictio = {};
  for listelt in listahan:
    dictio[cosetLabel(listelt[1])] = {};
    for dlgidl in listelt[2]:
      dictio[cosetLabel(listelt[1])][cosetLabel(dlgidl[1])] = convertIdeal(pariKr, sageKr, pariprelt_to_sage, dlgidl[2]);
  return(dictio)

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

# complexPolyToKrPoly, appears in 3.2
def complexPolytoKrPoly(Kr, pol, pr):
    CCpr = ComplexField(pr)
    emb = Kr.embedding(CCpr)
    newpol = recognize_polynomial(pol, Kr, emb=emb)
    return newpol

# naiveMatrixEqualityCheck, used by shimuraTheCosets
def naiveMatrixEqualityCheck(M1, M2, prec, tol):
  CCC = ComplexField(prec);
  MMM = MatrixSpace(CCC, 2);
  M1 = MMM(M1);
  M2 = MMM(M2);
  tol = tol;
  tentothenegtolerance = 10^-tol;
  matdif = M1-M2;
  dift1 = matdif[0,0];
  dift3 = matdif[1,0];
  dift2 = matdif[1,1];
  if abs(dift1.imag_part()) > tentothenegtolerance:
    return([None, 0]);
  if abs(dift2.imag_part()) > tentothenegtolerance:
    return([None, 0]);
  if abs(dift2.imag_part()) > tentothenegtolerance:
    return([None, 0]);
  if (abs(dift1.real_part() - round(dift1.real_part())) < tentothenegtolerance) & (abs(dift3.real_part() - round(dift3.real_part())) < tentothenegtolerance) & (abs(dift2.real_part() - round(dift2.real_part())) < tentothenegtolerance):
    return ([matdif, 1]);
  matsum = M1+M2;
  if (abs(CCC(matsum[0,0]).real_part() - round(CCC(matsum[0,0]).real_part())) < tentothenegtolerance) & (abs(CCC(matsum[1,0]).real_part() - round(CCC(matsum[1,0]).real_part())) < tentothenegtolerance) & (abs(CCC(matsum[1,1]).real_part() - round(CCC(matsum[1,1]).real_part())) < tentothenegtolerance):
    return ([matsum, -1]);
  return([None, 0]);

def translationMatrix(T, n=8):
  b1 = T[0,0].real_part().round();
  b3 = T[0,1].real_part().round();
  b2 = T[1,1].real_part().round();
  mm = Matrix( [ [1,0,b1,b3], [0,1,b3,b2], [0,0,1,0], [0,0,0,1] ] );
  uu = GSp_element(mat_convert(mm^-1, Zmod(8)));
  return(uu);

# shimuraTheCosets, appears in 2.2
# - cosets is an input
def shimuraTheCosets(cosets, Zs, bp, prec, tol):
  ret = {};
  Us = {};
  for coskey in cosets:
    cos = cosets[coskey];
    previ = None;
    for idlkey in cos:
      print(idlkey);
      idl = cos[idlkey];
      (U, m, u) = bp.Shimura_reciprocity(idl, n=8, period_matrix=True)
      origU = U;
      origu = u;
      yes = False;

      starti = previ;
      if starti == None:
        starti = 0;
      for i in range(starti, len(Zs)): # can be improved...
        Z = Zs[i];
        if yes:
          break;
        echec = naiveMatrixEqualityCheck(U.complex_matrix(), Z.complex_matrix(), prec, tol);
        if echec[1] != 0:
          yes = True;
          U = echec[1]*(i+1);
          if coskey not in ret.keys():
            ret[coskey] = {};
          u = u*translationMatrix(echec[0]);
          if previ == None:
            previ = i
          ret[coskey][idlkey] = [U, u, origU, origu]
      
      if yes == False:
        if coskey not in ret.keys():
          ret[coskey] = {};
        ret[coskey][idlkey] = [0, 0, origU, origu];
        print("          --> FAILED");
        # return(False);
  return(ret);

# cleverlyEvaluateAtInvariants, appears in 2.3
#   uses data from shimuradcosets to figure out which theta functions to evaluate
#   and then evaluates them
def cleverlyEvaluateAtInvariants(shimuradcosets, Zs, invs, prec, verbose=False):
  ret = {};
  for coskey in shimuradcosets:
    cos = shimuradcosets[coskey];
    if coskey not in ret.keys():
      ret[coskey] = {};
    for idlkey in cos:
      Uu = cos[idlkey];
      [Uind, u, origU, origu] = Uu;
      if Uind == 0:
        U = origU;
        u = origu;
      else:
        U = Zs[abs(Uind)-1];
      if idlkey not in ret[coskey].keys():
        ret[coskey][idlkey] = [0 for i in range(0, len(invs))];
      for i in range(0, len(invs)):
        if i != 0:
          continue;
        inv = invs[i];
        Unew = U;
        unew = u;
        print("Coset %5s  |  Ideal %5s  |  PM: Zs[%5d] |" % (coskey, idlkey, abs(Uind)-1), end="" );
        time_tmp = walltime(); # timing

        # Unew = Unew.complex_matrix(pr); #newline
        # choiceA = labrandeCalculThetaPari(inv^unew, Unew, pr); #newline
        choiceA = (inv^unew)(Unew,prec=prec,use_magma=True);
        choiceA2 = 1;

        if Uind < 0:
          choiceA = conjugate(choiceA);
          choiceA2 = conjugate(choiceA2);
        print(" time: %15.10f | %15.5f + %15.5fi | %15.5f + %15.5fi" % (walltime() - time_tmp, real_part(choiceA), imag_part(choiceA), real_part(choiceA2), imag_part(choiceA2) ) );
        ret[coskey][idlkey][i] = choiceA;
  return(ret);

# sanityCheck_smoothden
#   sets the GP variable F to the factorization of
#    the denominator of KrPoly
#   returns F 
def sanityCheck_smoothden(KrPoly):
  gp.set("F", gp.factor(KrPoly.denominator(), 10^6));
  return( [gp("F"), (gp("F[,1]~[#F[,1]]"))] );

# assumes smoothden is done
def sanityCheck_probsgal(pariKr, pariKrPoly, ub=100):
  return( (gp.isprobablygalois(pariKr, pariKrPoly, gp("Set(F[,1])"), ub) ) );

def sagePolyToPariPoly(poly, sageprelt_to_pari):
  polstr = str(poly).replace("y","x").replace("alphar","("+str(sageprelt_to_pari)+")")
  pol = gp(polstr)
  return(pol)

def rnfConductor(poly, pariKr, tulong=10000):
  cond = gp.rnfconductor(pariKr, [poly, tulong])
  if cond == 0:
    raise ValueError;
    return(0)
  return(cond[1])

# bleh = Theta_element_polynomial_ring object
# PM   = period matrix, typically Unew
# prec = precision

# PM = Zs[1].complex_matrix(2000);
# CCC = ComplexField(2000);
# z = Matrix([[CCC(0)],[CCC(0)]]);
# c = Matrix([[0],[0],[0],[0]]);

# zmag = magma(z)
# g = z.nrows()
# taumag = magma(PM)
# cmag = magma(c)
# magma.Theta(cmag, zmag, taumag).sage()

# magma.CalculThetas(zmag, taumag).sage()

# magma.load("fastthetasgenus2.m");
# magma.load("wamelentheta.m");
# bleh = Theta_element_polynomial_ring object
# PM   = period matrix, typically Unew
# prec = precision

# PM = Zs[1].complex_matrix(2000);
# CCC = ComplexField(2000);
# z = Matrix([[CCC(0)],[CCC(0)]]);
# c = Matrix([[0],[0],[0],[0]]);

# zmag = magma(z)
# g = z.nrows()
# taumag = magma(PM)
# cmag = magma(c)
# magma.Theta(cmag, zmag, taumag).sage()

# magma.CalculThetas(zmag, taumag).sage()

# magma.load("fastthetasgenus2.m");
# magma.load("wamelentheta.m");

@cached_function
def CalculThetas(zmag, taumag):
  return magma.CalculThetas(zmag, taumag, nvals=2);

@cached_function
def thetafunsq(i, zm, tm, pw):
  [fn, cn] = CalculThetas(zm, tm/2);
  return magma.thetafunsqhelper(converti(i), fn, cn).sage()^(pw/2);

def labrandeCalculTheta( bleh, PM, prec ):
    
    CCC    = ComplexField(ceil(prec*1.1));
    zmag   = magma(Matrix([[CCC(0)],[CCC(0)]]));
    taumag = magma(PM);

    numpol = bleh._num_pol.degrees();
    numcoe = bleh._num_pol.coefficients()[0];

    denpol = bleh._den_pol.degrees();
    dencoe = bleh._den_pol.coefficients()[0];

    retcoe = numcoe/dencoe;
    retnum = 1;
    retden = 1;

    for i in range(0, 16):
        if numpol[i]%2 != 0:
            error("SAD");
        if numpol[i] != 0:
            retnum *= thetafunsq(i, zmag, taumag, numpol[i]);
            # print("BOOM ");
        if denpol[i]%2 != 0:
            error("SAD");
        if denpol[i] != 0:
            retden *= thetafunsq(i, zmag, taumag, denpol[i]);
            # print("BAM ");

    return(retcoe*retnum/retden);


gp.read("fasttheta.gp");
gp.read("theta.gp");

@cached_function
def CalculThetasgp(zmag, taumag):
  print(taumag);
  return gp.CalculThetas(zmag, taumag);

@cached_function
def thetafunsqgp(i, zm, tm, pw):
  [fn, cn] = CalculThetasgp(zm, tm/2);
  return gp.thetafunsqhelper(converti(i), fn, cn).sage()^(pw/2);

def labrandeCalculThetaPari( fctrs, PM, prec ):
    
    CCC   = ComplexField(ceil(prec*1.1));
    zgp   = gp.mattranspose(gp([CCC(0),CCC(0)]));
    taugp = gp(PM);

    numpol = fctrs._num_pol.degrees();
    numcoe = fctrs._num_pol.coefficients()[0];

    denpol = fctrs._den_pol.degrees();
    dencoe = fctrs._den_pol.coefficients()[0];

    retcoe = numcoe/dencoe;
    retnum = 1;
    retden = 1;

    for i in range(0, 16):
        if numpol[i]%2 != 0:
            error("SAD");
        if numpol[i] != 0:
            retnum *= thetafunsqgp(i, zgp, taugp, numpol[i]);
            # print("BOOM ");
        if denpol[i]%2 != 0:
            error("SAD");
        if denpol[i] != 0:
            retden *= thetafunsqgp(i, zgp, taugp, denpol[i]);
            # print("BAM ");

    return(retcoe*retnum/retden);


def converti(i):
  # return i;
  arr = [0, 2, 1, 3, 8, 10, 9, 11, 4, 6, 5, 7, 12, 14, 13, 15];
  return arr[i];
# labrandeCalculTheta(inv^unew, Unew, prec)