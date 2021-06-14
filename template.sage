gp.set("SHIMURADIR", "\""+SHIMURADIR+"\"");
gp.read(FGAGDIR+"/fgag.gp");
gp.read(FGAGDIR+"/fgagshimuray.gp");
load("cm2functions.sage");
gp.read("theta.gp");

# Set precision.f

# m is the modulus
#   used to compute ShGrp
m = 2

# pr is the precision
#   used to compute theta functions
pr = 900*2*2*2#*2*2

gp.default("realprecision", floor(pr/4)-10);

# pr2 / tol
pr2 = 3000
tol = 10 # tolerance, number of digits, lower is more accurate but slower...

# Set DAB.
V = [17, 5, 2]
#V = [37, 13, 33]
# V = [41, 7, 2]
# V = [113, 33, 18]
# V = [17, 17, 34]
# V = [809, 53, 500]
# V = [28, 12, 29]

# V = [796, 52, 477]
# V = [5233, 75, 98]

# 0 Preliminaries

# 0.1 PariSetup
# - ShGrp, Shimura group as a PARI object
# - pariKr, bnfinit of Kr as a PARI object
# * first time m is used
time_tmp = walltime(); # timing
[ShGrp, pariKr] = parisetup(V, m=m)
print("%50s %10.10f" % ("PariSetup", walltime()-time_tmp) ) # timing

# 0.2 SageSetup
# - K, CM field
# - Zs, list of period matrices
# - bp = Zs[bpind]
# - sageKr, Kr as a SAGE object
# - yr, primitive element of sageKr
time_tmp = walltime(); # timing
[K, Zs, [bp, bpind], [sageKr, yr]] = sagesetup(V, 1);
print("%50s %10.10f" % ("SageSetup", walltime()-time_tmp) ); # timing

# 0.3 Primitive Elements
# - pariprelt, the primitive element of pariKr (i.e. pariKr = Kr(pariprelt))
# - sageprelt, the primitive element of sageKr (i.e. sageKr = Kr(sageprelt))
# pariKr and sageKr are related by an isomorphism
#   pariKr     -->  sageKr
#   pariprelt  |->  pariprelt_to_sage
# whose inverse is given by
#   sageKr     --> pariKr
#   sageprelt  |->  sageprelt_to_pari
# * first time pariKr is used
# * first time sageKr is used
time_tmp = walltime(); # timing
pariprelt = gp("y");
sageprelt = yr;
[pariprelt_to_sage, sageprelt_to_pari] = sagepariPrimitiveElement(pariKr, pariprelt, sageKr, sageprelt);
print("%50s %10.10f" % ("PrimElts", walltime()-time_tmp) ) # timing

# 1 Find G/H

# 1.1 PariCosets
# Let G = elements of class group which is isomorphic to Gal(CMKr(m)/Kr).
# Let H = elements of class group which is isomorphic to Gal(CMKr(m)/HKr(1)).
# This uses PARI to compute G/H = {gH : g in G} = Gal(HKr(1)/Kr),
# * first time ShGrp is used
time_tmp = walltime(); # timing
paricosets = pariCosets(ShGrp, m);
print("%50s %10.10f" % ("PariCosets", walltime()-time_tmp) ) # timing

# 1.2 SageCosets
# Convert paricosets to sagecosets, a dictionary.
# * first time pariprelt_to_sage is used
time_tmp = walltime(); # timing
sagecosets = sagepariIdealConversion(pariKr, sageKr, pariprelt_to_sage, paricosets)
print("%50s %10.10f" % ("SageCosets", walltime()-time_tmp) ) # timing

# 2 Theta Function Computations

# 2.1 Declare Rosenhain invariants
# Constant operation.
# May be changed up later.
time_tmp = walltime(); # timing
ros = rosenhain_invariants(2)
print("%50s %10.10f" % ("RosDeclarations", walltime()-time_tmp) ) # timing

# 2.2 Use Shimura reciprocity.
# Are pr2 and tol necessary?
time_tmp = walltime(); # timing
shimurad = shimuraTheCosets(sagecosets, Zs, bp, pr2, tol)
print("%50s %10.10f" % ("ShimuraTheCosets", walltime()-time_tmp) ) # timing

# 2.3 Cleverly evaluate at invariants
# Use the cosets to find the roots of a defining polynomial for HK(1).
# * uses shimurad, results of 2.2 for the first time
invs=rosenhain_invariants(2)
time_tmp = walltime(); # timing
dic = cleverlyEvaluateAtInvariants(shimurad, Zs, invs, pr, verbose=verbose)
print("%50s %10.10f" % ("CleverInvariantEvaluation", walltime()-time_tmp) ) # timing

# 3 The Polynomial

# 3.1 The Approximation
# From the theta function evaluations,
#   find a complex approximation of the polynomial we want.
time_tmp = walltime(); # timing
roots = [ sum([ dic[coskey][idlkey][0]  for idlkey in dic[coskey]]) for coskey in dic ]
complexPoly = rootstoPoly(roots)
print("%50s %10.10f" % ("ComplexRootsToPoly", walltime()-time_tmp) ) # timing

allroots = [ dic[key1][key2][0] for key1 in dic for key2 in dic[key1] ]
CM2polycomplex = rootstoPoly(allroots)

# 3.2 The Recognition
# Recognize the polynomial coefficients as an element of sageKr.
time_tmp = walltime(); # timing
CM2Poly = complexPolytoKrPoly(sageKr, CM2polycomplex, CM2polycomplex.coefficients()[0].prec())
print("%50s %10.10f" % ("Recognize", walltime()-time_tmp) ) # timing

# X = PolynomialRing(sageKr,'X');
CM2.<b> = sageKr.extension(CM2Poly)

# 4 The Checks

# 4.1 Factor the denominator.
#time_tmp = walltime(); # timing
#smallprime = 30000;
#while sanityCheck_smoothden(KrPoly)[1] > smallprime:
#	raise ValueError("INCREASE PRECISION? smooth denominator check failed: %s" % sanityCheck_smoothden(KrPoly)[0]);

# 4.2 Verify that it is not not Galois.
#pariKrPoly = sagePolyToPariPoly( KrPoly, sageprelt_to_pari )
#sanchk = gp.type( sanityCheck_probsgal(pariKr, pariKrPoly, ub=100) )
#if sanchk == "t_STR":
#	print(sanchk);
#	raise ValueError(sanchk);
#print("%50s %10.10f" % ("Sanity Checks", walltime()-time_tmp) ) # timing

# 4.3 Take a long time to actually compute the real conductor.
#time_tmp = walltime(); # timing
#cond = rnfConductor(pariKrPoly, pariKr, 10000)
#print("%50s %10.10f" % ("Conductor Correctness Check", walltime()-time_tmp) ) # timing
