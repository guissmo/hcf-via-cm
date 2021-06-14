/* NEW */

/* PRIVATE */

/* GENERIC GROUP */
xxxfgagmulgngp() = { return( ( (c, x, y)->(x+y) ) ); }
xxxfgagpowgngp() = { return( ( (c, x, n)->(n*x) ) ); }
xxxfgagdlggngp() = { return( ( (c, x)->lift(x) ) ); }
xxfgaginitgngp(bektor) = {

  my(ops);
  ops = [
    xxxfgagmulgngp(),
    xxxfgagpowgngp(),
    xxxfgagdlggngp()
  ];

  my(G, DG);
  G = vector(#bektor, i,
    vector(#bektor, j,
      if(bektor[j] == 0,
        i==j,
        Mod( i==j, bektor[j])
      )
    )~
  );
  DG = bektor;
  my(ret);
  ret = [ G, DG, ops, 0, 0];
  return(ret);
}

/* IDEAL/RAY CLASS GROUP */
xxxfgagmulclgp() = { return( ( (c, x, y)->idealmul(c, x, y)) ); }
xxxfgagpowclgp() = { return( ( (c, x, n)->if(n==0,return(idealhnf(c, 1)));idealpow(c, x, n)) ); }
xxxfgagdlgclgp() = { return( ( (c, x)->bnrispr=bnrisprincipal(c, x)[1]; if(#bnrispr == 0, return([1]~), return(bnrispr))) ); }
xxfgaginitclgp(bnr) = {

  if(#bnr==10, bnr=bnrinit(bnr, 1, 1));

  my(ops);
  ops = [
    xxxfgagmulclgp(),
    xxxfgagpowclgp(),
    xxxfgagdlgclgp()
  ];

  my(K, clgp);
  K = bnr;
  clgp = K.clgp;
  if(#K.gen == 0,
    clgp[3] = [idealhnf(K, 1)];
    clgp[2] = [1];
  );
  
  my(G, DG);
  G = clgp[3];
  DG = clgp[2];

  my(ret);
  ret = [ G, DG, ops, 1, bnr];
  return(ret);
}

/* UNIT GROUP */
xxxfgagmulungp() = { return( ((c, x, y)->nfeltmul(c, x, y))  ); }
xxxfgagpowungp() = { return( ((c, x, n)->if(n==0,return(1));nfeltpow(c, x, n))  ); }
xxxfgagdlgungp() = { return( ((c, x)->tmp=bnfisunit(c, x);tmp=concat([tmp[#tmp], tmp[1..(#tmp-1)]]);lift(tmp))  ); }
xxfgaginitungp(K:bnf) = {

  my(ops);
  ops = [
    xxxfgagmulungp(),
    xxxfgagpowungp(),
    xxxfgagdlgungp()
  ];

  my(G, DG);
  G = concat([K.tu[2], K.fu]);
  DG = concat([K.tu[1], vector(#K.fu)]);
  if(matsize(G) == [0, 0],
    G = Mat(1);
    DG = matdiagonal([1]);
  );

  my(ret);
  ret = [G, DG, ops, 2, K];
  return(ret);
}

/* IDEAL STAR */
xxxfgagmulidst() = { return( ((c, x, y)->nfeltmul(c[1], x, y))  ); }
xxxfgagpowidst() = { return( ((c, x, n)->if(n==0,return(1));nfeltpow(c[1], x, n))  ); }
xxxfgagdlgidst() = { return( ((c, x)->lift(ideallog(c[1], x, c[2])))  ); }
xxfgaginitidst(bnfm) = {

  my(bnf, m);
  [bnf, m] = bnfm;

  my(ops);
  ops = [
    xxxfgagmulidst(),
    xxxfgagpowidst(),
    xxxfgagdlgidst()
  ];

  my(K);
  K = bnf;

  my(bid);
  bid = idealstar(K, m, 2);
  
  my(G, DG);
  G = apply( x->nfbasistoalg(K, x), bid.gen);
  DG = bid.cyc;
  if(matsize(G) == [0, 0],
    G = Mat(1);
    DG = matdiagonal([1]);
  );

  my(ret);
  ret = [G, DG, ops, 3, [K, bid]];

  return(ret);
}

xxfgaggen(Gr) = { return(Gr[1]) };
xxfgagcyc(Gr) = { return(Gr[2]) };
xxfgagops(Gr) = { return(Gr[3]) };
xxfgagflg(Gr) = { return(Gr[4]) };
xxfgagcon(Gr) = { return(Gr[5]) };

/*************************/

fgaggen(Gr)       = { return( xxfgaggen(Gr) ) };
fgagcyc(Gr)       = { return( xxfgagcyc(Gr) ) };
fgagcard(Gr)      = { return( my(ret=vecprod(fgagcyc(Gr))); if(ret==0, return(+oo)); return(ret) ) };
fgagid(Gr)        = { return( fgagpow(Gr, 0, 0) ) };
fgagmul(Gr, x, y) = { return( (xxfgagops(Gr)[1])(xxfgagcon(Gr), x, y) ) };
fgagpow(Gr, x, n) = { return( (xxfgagops(Gr)[2])(xxfgagcon(Gr), x, n) ) };

fgagdlg(Gr, x, flag)    = {
  my(ret);
  ret = (xxfgagops(Gr)[3])(xxfgagcon(Gr), x);
  if(flag == 0 , return( ret ) );
  for(i = 1, #ret,
    my(ord);
    ord = fgagcyc(Gr)[i];
    if(ord != 0, ret[i] = Mod(ret[i], ord));
  );
  return(ret);
};

fgaginit(context, flag) = {
  if(flag == 0, return( xxfgaginitgngp(context) )); /* generic matrix group */
  if(flag == 1, return( xxfgaginitclgp(context) )); /* class group of bnf / bnr */
  if(flag == 2, return( xxfgaginitungp(context) )); /* unit group */
  if(flag == 3, return( xxfgaginitidst(context) )); /* ideal star */
}

/**************************/

fgagsubgpintsum(Gr, H1, H2, intsumflag) = {

  my(B, DB);
  [B, DB] = [fgaggen(Gr), matdiagonal(fgagcyc(Gr))];

  my(H1H2, H12, U12);
  H1H2 = matconcat([H1, H2]);
  [H12, U12] = mathnf(H1H2, 1);

  my(U1, U3);
  U1 = U12[1..#B, 1..#B];
  U3 = U12[(#B+1)..#U12, 1..#B];

  my(H3, H);
  H3 = mathnf(matconcat([H1*U1, DB])); /* intersection */
  H = mathnf(matconcat([H2*U3, DB]));  /* sum */

  if(intsumflag == 1, return(H3));
  if(intsumflag == 2, return(H));
  return([H3, H]);

}

fgagsubgpint(Gr, H1, H2) = { return(fgagsubgpintsum(Gr, H1, H2, 1)); }
fgagsubgpsum(Gr, H1, H2) = { return(fgagsubgpintsum(Gr, H1, H2, 2)); }

/* checks if H1 is subgroup of H2 */
fgagissubgp(Gr, H1, H2) = {
  my(Hint);
  Hint = fgagsubgpint(Gr, H1, H2);
  if(Hint == H1, return(1));
  return(0);
}

/**************************/

xxxfgagsnfdlg(dlg, Ua, DG) = {
  return ( (c,x)->
                  d=Ua*dlg(c,x);
                  for(i = 1, #d,
                    if(DG[i] != 0,
                      d[i] = lift(Mod(d[i], DG[i]));
                    );
                  );
                  return(d);
  )
};

xxfgagmulbasis(Gr, B, M) = {

  my(nM, mM, ret);
  [nM, mM] = matsize(M);
  if( #B != nM, error("dimension error", str) );

  ret = vector(mM, i, fgagid(Gr) );

  for(j = 1, mM,
    for(k = 1, nM,
      ret[j] = fgagmul( Gr, ret[j], fgagpow( Gr, B[k], M[k, j] ) );
    );
  );
  return(ret);

}

xxfgagsubgprels(HG, Gr) = {

  my(DG);
  DG = matdiagonal(fgagcyc(Gr));

  my(GHG);
  GHG = xxfgagmulbasis( Gr, fgaggen(Gr), HG );

  my(almostSNF);
  almostSNF = [GHG, HG^-1*DG];

  my(newdlg);
  newdlg = (xxxfgagsnfdlg(Gr[3][3], HG^-1, fgagcyc(Gr)));

  return( [GHG, HG^-1*DG, [Gr[3][1], Gr[3][2], newdlg], -1, xxfgagcon(Gr)] );

}

fgagsubgpsnf( HG, Gr ) = {

  if(fgagcard(Gr) == +oo, error("group is not finite") );

  my(GM);
  GM = xxfgagsubgprels(HG, Gr); /* not yet tested if works ... apparently it does */

  my(SNF);
  SNF = fgagsnf( GM[1..2], GM );

  return(SNF);

}

xxfgagtrivopercxy( id ) = {
  return( (c,x, y)->id );
}

xxfgagtrivopercx( id ) = {
  return( (c,x)->id );
}

xxfgaginittrivgrp( id ) = {
  \\ warning("trivial group ahead!");
  return(
  [
    [id],
    [1],
    [
      xxfgagtrivopercxy(id),
      xxfgagtrivopercxy(id),
      xxfgagtrivopercx([1]~)
    ],
    -1,
    -1
  ]
  );
}

fgagsnf( GM, Gr ) = {

  my(G, M);
  [G, M] = GM;

  my(H);
  H = mathnf(M);

  my(nH, mH);
  [nH, mH] = matsize(H);
  if( nH != mH, error("H is not a square matrix / M is not of maximal rank OR the group is infinite", str) );

  my(U, V, D);
  [U, V, D] = matsnf(H, 1);

  my(nD, mD);
  [nD, mD] = matsize(D);
  if( nD != mD, error("D is not a square matrix", str) );

  my(Ap);
  Ap = xxfgagmulbasis(Gr, G, U^-1);

  my(n);
  for(i=1, nD,
    if(D[i,i] == 1, break());
    n = i;
  );

  my(A, DA, Ua);
  A = Ap[1..n];
  DA = vector(n, i, D[i, i]);
  Ua = U[1..n,];
  if(vecprod(matsize(A)) == 0, return( xxfgaginittrivgrp(fgagid(Gr)) ));

  my(mul, pow, dlg, mulpowdlg, newdlg);
  [mul, pow, dlg] = xxfgagops(Gr);
  newdlg = xxxfgagsnfdlg(dlg, Ua, DA);

  /* if(#A == 0, A = [fgagid(Gr)]; DA = [1] ); */

  my(ret);

  ret = [A, DA, [mul, pow, newdlg], -1, xxfgagcon(Gr)];
  return(ret);

}

/**************************/

/* SECTION: MORPHISMS AS INPUT */

/* Quotient group B/A of two (normal) subgroups A and B of a group G */
/* Assuming A is a subgrp of B */
/* Thus psi = identity */
fgagsubgpquo(Gr, HA, HB) = {

  my(GrQuo, BH, HinvDB, mulpowdlg, mul, pow, trash, dlg, context);
  GrQuo = [BH, HinvDB, mulpowdlg, trash, context] = xxfgagsubgprels(HB, Gr);
  mulpowdlg = [mul, pow, dlg];

  my(G, M);
  G = BH;
  M = HA*HB^-1;

  my(SNF);
  SNF = fgagsnf( [G, M], GrQuo);

  return( SNF );

}

/**************************/

/* inverse image */
fgagmorsubgpinvimg(GrB, GrC, morBtoC, HC, dlgflag) = {

  my(B, DB); my(C, DC);
  [B, DB] = [fgaggen(GrB), matdiagonal(fgagcyc(GrB))];
  [C, DC] = [fgaggen(GrC), matdiagonal(fgagcyc(GrC))];

  my(morBtoCdlog);
  morBtoCdlog = morBtoC;
  if( dlgflag == 0, morBtoCdlog = ((x)->fgagdlg(GrC, morBtoC(x))) );

  my(P);
  P = matconcat(apply(morBtoCdlog, B));

  my(HC);
  if(HC == 0, HC = DC);

  my(nP, mP, nHC, mHC);
  [nP, mP] = matsize(P);
  [nHC, mHC] = matsize(HC);
  /* sanity check */
  if(nP != nHC, error("incompatible size for P and HC", str) );

  my(PHC);
  PHC = matconcat([P, HC]);

  my(H, U);
  [H, U] = mathnf(PHC, 1);

  /* U1 is a square matrix with dimensions nP x mP */
  my(U1);
  U1 = U[1..mP,1..mP];

  my(U1DB, HB);
  U1DB = matconcat( [U1, DB] );
  HB = mathnf(U1DB);

  return(HB);

  if(hnfflag == 0, return(HB));
  return( [ HB, GrB ] );

}

fgagmorsubgpimg(GrB, GrC, morBtoC, HB, dlgflag) = {

  my(B, DB); my(C, DC);
  [B, DB] = [fgaggen(GrB), matdiagonal(fgagcyc(GrB))];
  [C, DC] = [fgaggen(GrC), matdiagonal(fgagcyc(GrC))];

  my(morBtoCdlog);
  morBtoCdlog = morBtoC;
  if( dlgflag == 0, morBtoCdlog = ((x)->fgagdlg(GrC, morBtoC(x))) );

  my(P);
  P = matconcat(apply(morBtoCdlog, B));
  
  my(PHB, M);
  PHB = P*HB;
  M = matconcat([PHB, DC]);

  my(HC);
  HC = mathnf(M);
  return(HC);

  if(hnfflag == 0, return(HC));
  return( [ HC, GrC ] );

}

fgagmorker(GrB, GrC, morBtoC, dlgflag) = { return( fgagmorsubgpinvimg(GrB, GrC, morBtoC, fgagcyc(GrC), dlgflag)); }
fgagmorimg(GrB, GrC, morBtoC, dlgflag) = { return( fgagmorsubgpimg(GrB, GrC, morBtoC, matid(#fgaggen(GrB)), dlgflag) ); }

/**************************/

/* assuming c is context */
xxgrpextGrA(c) = { return(c[1]) };
xxgrpextGrC(c) = { return(c[2]) };
xxgrpexteff(c) = { return(c[3]) };
xxgrpexteffinv(c) = { return(c[4]) };
xxgrpextgee(c) = { return(c[5]) };
xxgrpextgeeinv(c) = { return(c[6]) };
xxgrpextinputcon(c) = { return(c[7]) };
xxgrpextBp(c) = { return(c[8]) };
xxgrpextGrBmul(c) = { return(c[9][1]) };
xxgrpextGrBpow(c) = { return(c[9][2]) };

/*
xxxfgaggrpextdlg() = {

  return(
    (
      (c, x)->

         my(id, mul, pow, red, gee, effinv, Bp);
         mul = agc_grpext_GrB_mul(c);
         pow = agc_grpext_GrB_pow(c);
         gee = agc_grpext_GrB_gee(c);
         effinv = agc_grpext_GrB_effinv(c);
         Bp = agc_grpext_GrB_Bp(c);
         
         my(Y);
         Y = ag_dlg(agc_grpext_GrC(c), gee(c, x));

         my(BpY); /* BpY is an element of B.
         BpY = ag_mulbasis([0, 0, [mul, pow, 0], -1, c], Bp, Mat(Y));
         if(matsize(BpY) != [1,1], error("error"));
         BpY = BpY[1,1];

         my(xdivBpY);
         xdivBpY = mul(c, x, pow(c, BpY, -1));

         my(X);
         X = ag_dlg(agc_grpext_GrA(c), effinv(c, xdivBpY));

         return(concat([X~,Y~])~);
    )
  );
}
*/

xxxfgaggrpextdlg() = {

  return(
    (
      (c, x)->
      my(oldcon);
      oldc = xxgrpextinputcon(c);
      my(Y);
      Y = fgagdlg(xxgrpextGrC(oldc), xxgrpextgee(c)(oldc, x));

      my(BpY); /* BpY is an element of B. */
      BpY = xxfgagmulbasis([0, 0, [xxgrpextGrBmul(c), xxgrpextGrBpow(c), 0], -1, oldc], xxgrpextBp(c), Mat(Y));
      if(#BpY != 1, error("error"));
      BpY = BpY[1];

      my(xdivBpY);
      xdivBpY = xxgrpextGrBmul(c)(oldc, x, xxgrpextGrBpow(c)(oldc, BpY, -1));

      my(X);
      X = fgagdlg(xxgrpextGrA(oldc), xxgrpexteffinv(c)(oldc, xdivBpY));

      return(concat([X~,Y~])~);
    )
  );
}


/*
  This function returns an fgag structure
    representing the middle term of a short exact sequence
      of the form
        0 --> A --f--> B --g--> C --> 0
    given
      GrA       : fgag structure representing GrA
      GrC       : fgag structure representing GrC
      eff       : a closure of the form (c,x)->dostuff
                    which represents the function f: A -> B
      effinv    : a closure of the form (c,y)->dostuff
                    which represents a section of f
                      that is, it outputs x such that f(x) = y
      gee       : a closure of the form (c,x)->dostuff
                    which represents the function g: B -> C
      geeinv    : a closure of the form (c,y)->dostuff
                    which represents a section of g
                      that is, it outputs x such that g(x) = y
      mulpowcon : a vector of length three where
                    mulpowcon[1]
                      is a closure of the form (c,x,y)
                        which represents the mul operation in B
                    mulpowcon[2]
                      is a closure of the form (c,x,n)
                        which represents the pow operation in B
                    mulpowcon[3]
                      is a vector
                        which contains context for the following
                        closures which involve B:
                          eff
                          effinv
                          gee
                          geeinv
                          mul
                          pow
*/
fgaggrpext(GrA, GrC, eff, effinv, gee, geeinv, mulpowcon, flag) = {

  my(A, DA, C, DC);
  [A, DA] = [fgaggen(GrA), matdiagonal(fgagcyc(GrA))];
  [C, DC] = [fgaggen(GrC), matdiagonal(fgagcyc(GrC))];

  my(mul, pow, c);
  [mul, pow, c] = mulpowcon;

  /* Find Bp such that gee(Bp) = C. */
  my(Bp);
  Bp = apply( x->geeinv(c,x), C);

  /* Set Bpp = Bp*DC. */
  my(Bpp);
  Bpp = xxfgagmulbasis([0, 0, [mul, pow, 0], -1, c], Bp, DC);

  /* Find App such that eff(App) = Bpp. */
  my(App);
  App = apply( x->effinv(c,x), Bpp);

  my(P);
  P = matconcat( apply( (A_elt)->fgagdlg(GrA, A_elt), App ) );
  
  if(flag == 1, return(P));

  /* Find effA. */
  my(effA);
  effA = apply( x->eff(c,x), A);

  my(G, M);
  /* G = [ f(A) | Bp ] */
  G = concat([effA, Bp]);
  M = matconcat([DA, -P; 0, DC]);
  /*breakpoint();*/

  my(dlgwrtGM);
  dlgwrtGM = xxxfgaggrpextdlg();

  if(flag == 2, return([G, M, [mul, pow, dlgwrtGM], 1, c]));

  my(SNF);
  SNF = fgagsnf([G, M], [0, 0, [mul, pow, dlgwrtGM], -1, c]);

  my(newcon);
  newcon = [GrA, GrC, eff, effinv, gee, geeinv, c, Bp, [mul, pow] ];

  /* update con */
  SNF[5] = newcon;

  my(retGr);
  retGr = SNF;

}

/**************************/

fgaglisteltnaive(Gr, flag) = {
  my(ret);
  ret = List();
  if(Gr[1] == [], listput(/*~*/ret, fgagid(Gr)); return(ret));
  forvec(X = vector(#fgagcyc(Gr), i, [0, fgagcyc(Gr)[i]-1]),
    my(elt);
    elt = xxfgagmulbasis(Gr, fgaggen(Gr), Mat(X~))[1];
    if(flag == 1,
      listput(/*~*/ret, [X, elt]),
      listput(/*~*/ret, elt)
    );
  );

  return(ret);
}

fgaglistcoseltnaive(Gr, ClKrm, coset, elt, flag) = {
  my(ret);
  if(ClKrm == 0, ClKrm = Gr);
  ret = coset;
  for(i = 1, #ret, 
    if(flag == 1,
      my(a = fgagmul(Gr, coset[i][2], elt));
      ret[i] = [fgagdlg(ClKrm, a, 1), a],
      ret[i] = fgagmul(Gr, coset[i], elt);
    )
  );
  return(ret);
}

fgaglistcoseltsnaive(Gr, ClKrm, H, flag) = {
  my(Quo);
  Quo = fgagsubgpquo(Gr, H, matid(#H));
  my(Hsnf, Hlist, GHreps);
  Hsnf = fgagsubgpsnf(H, Gr);
  Hlist = fgaglisteltnaive(Hsnf, flag);
  GHreps = fgaglisteltnaive(Quo, flag);
  my(ret);
  ret = List();
  for(i=1, #GHreps,
    if(flag == 1,
      my(dlooog);
      dlooog = vector( #fgagcyc(Quo), j, Mod( GHreps[i][1][j] , fgagcyc(Quo)[j] ) );
      listput(/*~*/ret, [ dlooog, fgaglistcoseltnaive(Gr, 0, Hlist, GHreps[i][2], 1) ]),
      listput(/*~*/ret, fgaglistcoseltnaive(Gr, 0, Hlist, GHreps[i], 1));
    )
  );
  return(ret);
}
