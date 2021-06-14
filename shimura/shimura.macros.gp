/*
listcmfields() = {
	v = externstr("ls -ldtr *.pol | gawk -F\" \" '{print $9}' | gawk -F. '{print $1}' | gawk -F_ '{print $1}{print $2}{print $3}'");
	v = apply(x->eval(x), v);
	return ( vector(#v/3, i, [  v[(i-1)*3 + 1], v[(i-1)*3 + 2], v[(i-1)*3 + 3]  ]) );
}
*/

cmfield(D, A, B) = {
	if(#D==3, [D,A,B]=D;return(init_cmfield(A,B)));
	if(#D==2, [A,B]=D;return(init_cmfield(A,B)));
	if(#D==1, return(init_cmfield(A,B)));
	return(0);
}

cm_K0 = cmfield_K0(cm) = { return(cm[1]) }
cm_K = cmfield_K(cm) = { return(cm[2]) }
cm_Krel = cmfield_Krel(cm) = { return(cm[3]) }
cm_F0 = cm_K0r = cmfield_F0 = cmfield_K0r(cm) = { return(cm[6]) }
cm_F = cm_Kr = cmfield_F = cmfield_Kr(cm) = { return(cm[7]) }
cm_Krrel = cmfield_Krrel(cm) = { return(cm[8]) }
cm_DAB = cmfield_DAB(cm) = {
	my(p = cmfield_K(cm).pol);
	A = polcoeff(p, 2);
	B = polcoeff(p, 0);
	D = coredisc(A^2 - 4*B);
	return([D,A,B]);
}
cm_DABr(cm) = {
	my(p = cmfield_K(cm).pol);
	A = cm[14];
	B = cm[15];
	D = coredisc(A^2 - 4*B);
	return([D,A,B]);	
}

cm_print(cm) = {
	my(cm2 = cmfield(cm_DABr(cm)));
	printf("=======================================================================\n");
	printf("      %30s    %30s\n",cm_DAB(cm),cm_DABr(cm));
	printf("-----------------------------------------------------------------------\n");
	printf("  SG: %30s    %30s\n",cm[16][1],cm2[16][1]);
	printf(" TNS: %30s    %30s\n",cm[17][1],cm2[17][1]);
	printf(" CG : %30s    %30s\n",cm_K(cm).clgp[2],cm_K(cm2).clgp[2]);
	printf(" CG0: %30s    %30s\n",cm_K0(cm).clgp[2],cm_K0(cm2).clgp[2]);
	printf("  h : %30s    %30s\n",cm_K(cm).clgp[1],cm_K(cm2).clgp[1]);
	printf("  h0: %30s    %30s\n",cm_K0(cm).clgp[1],cm_K0(cm2).clgp[1]);
	printf(" igu: %30ld    %30ld\n",vecprod(cm2[17][1]),vecprod(cm[17][1]));
	printf("=======================================================================\n");
}


/* cm_intermediate */

/* finds complex conjugate ideal */
idealbar(K, af) = {
	return(complex_conjugate_ideal(K, idealhnf(K, af)));
}

/* macros */
conjmatrix(cm) = {
	return(shimura_group_precompute_conjugation_action(cm));
}

symplecticbases(cm) = {
	return(symplectic_bases(cm));
}

/* Given an ideal from idealhnf, and the corresponding xi, this expresses the symplectic basis of this ideal. */
symplectic_basis_unconverted_ideal(K, id, xi) = {
	return(symplectic_basis(K, switch_ideal_basis_for_symplectic_input(K,id), xi));
}

/* Given an ideal of bnf, determine the "bnfisprincipal" of aabard*/
aabardiff(bnf, xxx, bnr) = {
    return(idealmul(bnf,bnf.diff,idealmul(bnf,complex_conjugate_ideal(bnf,xxx),xxx)));
    if(bnr==0,
      return(bnfisprincipal(bnf,idealmul(bnf,bnf.diff,idealmul(bnf,complex_conjugate_ideal(bnf,xxx),xxx)))[1]),
      return(bnrisprincipal(bnr,idealmul(bnf,bnf.diff,idealmul(bnf,complex_conjugate_ideal(bnf,xxx),xxx)))[1]),
    );
}

/* Only works for dihedral case because I'm too lazy and we don't need the Galois case! */
typenorm_Kr_to_K(cmfield, ar) = {
   my (Kr, Lrel, Lrrel, a, n);

   Kr = cmfield [7];
   L_over_K = cmfield [11];
   L_over_Kr = cmfield [12];

   a = rnfidealnormrel(L_over_K, rnfidealup (L_over_Kr, ar));
   n = idealnorm (Kr, ar);

   return ([a, n]);
};

/* Only works for dihedral case because I'm too lazy and we don't need the Galois case! */
typenorm_K_to_Kr(cmfield, ar) = {
   my (Kr, Lrel, Lrrel, a, n);

   Kr = cmfield [7];
   L_over_K = cmfield [11];
   L_over_Kr = cmfield [12];

   a = rnfidealnormrel(L_over_Kr, rnfidealup (L_over_K, ar));
   n = idealnorm (K, ar);

   return ([a, n]);
};

/*
  Returns a matrix whose columns are the symplectic basis elements.
  The symplectic basis elements are expressed with respect to the power basis.
*/
symplecticbasis_wrt_powerbasis(K, af, xi) = {
	return( symplectic_transformation (symplectic_form (K, af, xi)) );
}

/*
  Given a period matrix M whose elements are from Kr,
  Embeds the matrix in the correct embedding, and thus resulting to a
    period matrix whose entries are complex numbers.
*/
embedperiodmatrix(cmKr, M) = {
    my(Kr, ret);
    if(#cmKr > 11, Kr = cm_Kr(cmKr), Kr = cmKr;);
    ret = 0;
	for(i=1,2,
    for(j=1,2, /* signs */
  	  my(a, b, c, E);
      if(j==2, M = subst(M, y, -y));
  	  a = nfeltembed(Kr, M[1,1])[i];
  	  b = nfeltembed(Kr, M[1,2])[i];
  	  c = nfeltembed(Kr, M[2,2])[i];
  	  E = [a, b; b, c];
  	  iE = imag(E);
  	  if( iE[1,1] > 0,
  	    if( matdet(iE) >= 0,
  	      ret = E;
  	    )
  	  )
    )
	);
  if(ret == 0, print("sad life"));
	return( ret );
}

periodmatrixpath = "/home/guissmo/math/mynotes/code/data/";

writeperiodmatrix(cm, M, prec, ind, prefix) = {
    if(prec < 20, default(realprecision, 20), default(realprecision, prec));
    if(prefix == 0, prefix = "");
    my(mat, D, A, B, fname, file);
    mat = embedperiodmatrix(cm, M);
    [D, A, B] = cm_DAB(cm);
    fname = concat([periodmatrixpath, prefix, A, "_", B, "_", ind, ".pm"]);
    file = fileopen(fname, "w");
	filewrite1(file, Str(real(mat[1,1])) ); filewrite1(file, " " );
	filewrite(file, Str(imag(mat[1,1])) );
	filewrite1(file, Str(real(mat[1,2])) ); filewrite1(file, " " );
	filewrite(file, Str(imag(mat[1,2])) );
	filewrite1(file, Str(real(mat[2,2])) ); filewrite1(file, " " );
	filewrite(file, Str(imag(mat[2,2])) );
	fileclose(file);
}

/*
  Given a cm structure, a matrix M = [tau1, tau2; tau2; tau3] over Kr,
  This writes, line-by-line, tau_i/2's on the file halfA_B_ind.pm.
*/
writehalfperiodmatrix(cm, M, prec, ind) = {
  writeperiodmatrix(cm, M/2, prec, ind, "half");
}

/*
  ONLY WORKS WHEN Kr.clgp IS CYCLIC.
  Given a cm structure, finds all period matrices, divides by 2, embeds them appropriately
    and then writes them on files halfA_B_ind.pm.
*/
writehalfperiodmatrices(cm, prec) = {
  pm = period_matrices(cm);
  pm = pm[1];
  for(i=1, #pm, writehalfperiodmatrix(cm, pm[i][1], prec, i));
}

padvec(v, l) = {
	return( concat(v, vector(l - #v)) );
}

intbasis_wrtpolybasis(K) = {
	return (Mat(vector(4, i, padvec(Vecrev(Vec(K.zk[i])),4)~)));
}

theta10_to_theta16(thvec) = {
  my(ret = vector(16));
  my(thetalegend = [1, 2, 3, 4, 6, 8, 9, 12, 15]);
  for(i=1, #thetalegend,
    ret[thetalegend[i]] = thvec[i+1];
  );
  ret[16] = thvec[1];
  return(ret);
}

theta10(cm) = {
  writehalfperiodmatrices(cm, 100);
  [D, A, B] = cm_DAB(cm);
  thvec = extern(concat(["/home/guissmo/math/mynotes/code/execs/igusa_streng -t < /home/guissmo/math/mynotes/code/data/half", Str(A), "_", Str(B), "_", Str(1), ".pm"]));
  return(theta10_to_theta16(thvec));
}

theta16_from_tau(mat) = {
  mathalf = mat/2;
  fname = concat([periodmatrixpath, "tmp.pm"]);
  file = fileopen(fname, "w");
  filewrite1(file, Str(real(mathalf[1,1])) ); filewrite1(file, " " );
  filewrite(file, Str(imag(mathalf[1,1])) );
  filewrite1(file, Str(real(mathalf[1,2])) ); filewrite1(file, " " );
  filewrite(file, Str(imag(mathalf[1,2])) );
  filewrite1(file, Str(real(mathalf[2,2])) ); filewrite1(file, " " );
  filewrite(file, Str(imag(mathalf[2,2])) );
  fileclose(file);
  thvec = extern(concat(["/home/guissmo/math/mynotes/code/execs/igusa_streng -t < ", fname]));
  return(theta10_to_theta16(thvec));
}

/* Output the Rosenhain invariants given the values of the 16 theta functions at tau. */
rosenhain(tau) = {
  theta16 = theta16_from_tau(tau);
  tupa1 = theta16[16]*theta16[2]/theta16[3]/theta16[1];
  tupa2 = theta16[2]*theta16[12]/theta16[1]/theta16[15];
  tupa3 = theta16[16]*theta16[12]/theta16[3]/theta16[15];
  return( [tupa1, tupa2, tupa3] );
}

rosenhain_from_theta16(theta16) = {
  tupa1 = theta16[16]*theta16[2]/theta16[3]/theta16[1];
  tupa2 = theta16[2]*theta16[12]/theta16[1]/theta16[15];
  tupa3 = theta16[16]*theta16[12]/theta16[3]/theta16[15];
  return( [tupa1, tupa2, tupa3] );
}

tau_from_Krideal(cm, frakb_from_Kr) = {
  frakb = frakb_from_Kr;
  tn = typenorm_Kr_to_K(cm, frakb);
  nB = idealnorm(K, frakb);
  bp = find_basepoint(cm, y);
  sb = symplectic_basis_unconverted_ideal(K, tn[1], xibp/nB);
  pm = period_matrix(cm, sb);
  em = embedperiodmatrix(cm, pm);
  return(em);
}

basa(s) = {
	read(concat("~/math/mynotes/code/",s));
}

algbasis_to_matrix(M) = {
  return ( Mat(apply(x->padvec(Vecrev(x),4)~,M))~ );
}

/* NICE FUNCTIONS*/

/*
  Given a bnr, idl_input, bnf:
  - Finds an ideal 
*/
findxi(bnr, idl_input, bnf, flag) = {
  my(idl, idl_dlog, idl_gen);
  idl = idealhnf(bnf, idl_input);
  Daabarinv = idealinv(bnf, idealmul(bnf, bnf.diff, idealmul(bnf, idl, complex_conjugate_ideal(bnf, idl))));
  [idl_dlog, idl_gen] = bnrisprincipal(bnr, Daabarinv);
  if(idl_dlog != 0, return(0));
  my(xi);
  xi = bnf.zk*idl_gen;
  un = [1, -1, lift(K.fu[1]), -lift(K.fu[1])];
  for(u=1,#un,
    xiun = nfeltmul(K, xi, un[u]);
    embs = nfeltembed(Kr, xiun);
    embs_imag = apply(x->imag(x), embs);
    if(embs_imag[1] < 0, next());
    if(embs_imag[2] < 0, next());
    if(flag == 1, Mod(xiun, bnf.pol));
    return( idl_gen );
  );
  return(0);
}

