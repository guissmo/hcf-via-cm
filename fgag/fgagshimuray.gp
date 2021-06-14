/* fgaginit */
/*~/pari-git/gp -f fgag.gp*/
SHIMURADIR="/home/guissmo/projects/hcfviacm/shimura";
read(concat(SHIMURADIR,"/shimura.gp"));
read(concat(SHIMURADIR,"/shimura.macros.gp"));

cm_init_from_twopols(Kpol, Krpol) = {

  my(A, B, Ar, Br);
  A = polcoeff(Kpol, 2);
  B = polcoeff(Kpol, 0);
  Ar = polcoeff(Krpol, 2);
  Br = polcoeff(Krpol, 0);

  my(K0, KoverK0, K);
  K0 = bnfinit(z^2 + A*z + B, 1);
  KoverK0 = rnfinit(K0, y^2 - z);
  K = bnfinit(y^4 + A*y^2 + B, 1);

  my(K0r, KroverK0r, Kr);
  K0r = bnfinit(z^2 + Ar*z + Br, 1);
  KroverK0r = rnfinit(K0r, y^2 - z);
  Kr = bnfinit(y^4 + Ar*y^2 + Br, 1);

  my(LabspoloverK, LabspoloverKr, LoverK, LoverKr);
  LabspoloverK = nffactor( K, substpol(nfsplitting(K.pol), y, x) )[1, 1];
  LabspoloverKr = nffactor( Kr, substpol(nfsplitting(Kr.pol), y, x) )[1, 1];
  LoverK = rnfinit(K, LabspoloverK);
  LoverKr = rnfinit(Kr, LabspoloverKr);

  my(sigm);
  if( poldegree(LoverK.mypolabs) == 8, sigm = 0, sigm = nfgaloisconj(Kr)[3] );

  my(ret);
  ret = vector(15);
  
  ret[ 1] = K0;
  ret[ 2] = K;
  ret[ 3] = KoverK0;
  
  ret[ 6] = K0r;
  ret[ 7] = Kr;
  ret[ 8] = KroverK0r;

  ret[11] = LoverK;
  ret[12] = LoverKr;
  ret[13] = sigm;

  return(ret);

}

cm_from_whatever(tmp) = {
  
  my(A, B);
  [A, B] = tmp[1..2];

  if(#tmp == 2, 
    if( type(tmp[1]) == "t_INT",
      [A, B] = [tmp[1], tmp[2]], /* from [A, B] */
      return( cm_init_from_twopols(A, B) ); /* from [Kpol, Krpol] */
    )
  ); 

  [A, B] = [-1, -1];

  if(#tmp == 3, [A, B] = [tmp[2], tmp[3]]); /* from [D, A, B] */

  my(pol);
  if(#tmp == 9,
    pol = tmp.pol;
    [A, B] = [polcoeff(pol, 2), polcoeff(pol, 0)];
  ); /* from nfinit */
  if(#tmp == 10,
    pol = tmp.pol;
    [A, B] = [polcoeff(pol, 2), polcoeff(pol, 0)];
  ); /* from bnfinit */

  if(#tmp >= 15, return(tmp));

  if(A == -1, error("invalid input for cm_from_whatever"));
  if(B == -1, error("invalid input for cm_from_whatever"));

  my(cm);
  cm = init_cmfield_basic(A, B);
  return(cm);

}


embedsign(F, x) = {
  
  xemb = nfeltembed(F, x);
  return ( sign(real(xemb[1])) + sign(real(xemb[2])) );

}

/***/
shraycontextGrA(c) = { return(c[1]) };
shraycontextGrC(c) = { return(c[2]) };
shraycontextK(c) = { return(c[3]) };
shraycontextKoverK0(c) = { return(c[4]) };
shraycontextKm(c) = { return(c[5]) };
shraycontextK0plus(c) = { return(c[6]) };


shraymul() = { 
  return(
    ( (c,x,y)->
      if(#c==9, c = c[7]);
      return([
     	  idealmul(shraycontextKm(c), x[1], y[1]),
        nfeltmul(shraycontextK0plus(c), x[2], y[2])
      ])
    );
  );
}

shraypow() = { 
  return(
    ( (c,x,n)->
      if(#c==9, c = c[7]);
      if(n == 0, return( [idealhnf(shraycontextKm(c), 1), 1] ));
      return([
        idealpow(shraycontextKm(c), x[1], n),
        nfeltpow(shraycontextK0plus(c), x[2], n)
      ]);
    );
  );
}

shrayeff() = {
  return (
    ( (c, A_elt)->
      [
        idealhnf(
          shraycontextK(c),
          rnfeltup(shraycontextKoverK0(c), A_elt)
        ),
        A_elt
      ]
    );
  );
}

shraygee() = {
  return (
    ( (c, B_elt)->
      B_elt[1]
    )
  );
}

shrayeffinv() = {
  return (
    ( (c, B_elt)->
      
      my(elt_div_BpY_dlog, elt_div_BpY_gen);
      [elt_div_BpY_dlog, elt_div_BpY_gen] = bnrisprincipal(shraycontextKm(c), B_elt[1]);
      if( #elt_div_BpY_dlog == 0, elt_div_BpY_dlog = 0);

      my(x_of_xOK, xbarx);
      x_of_xOK = elt_div_BpY_gen; /* as element of K */
      xbarx = rnfeltnorm(shraycontextKoverK0(c), x_of_xOK);

      my(red_B_elt);
      red_B_elt = [
        fgagmul(
          shraycontextGrC(c),
          B_elt[1],
          fgagpow(
            shraycontextGrC(c),
            rnfeltup(shraycontextKoverK0(c), x_of_xOK),
            -1
          )
        ),
        fgagmul(
          shraycontextGrA(c),
          B_elt[2],
          fgagpow(
            shraycontextGrA(c),
            xbarx,
            -1
          )
        )
      ];

      my(redfrakb, redbeta);
      [redfrakb, redbeta] = red_B_elt;

      if( redfrakb != fgagid(shraycontextGrC(c)), error("eff^-1 is incompatible") );
      return(redbeta);
    )
  )
}

shraygeeinv() = {
  return (
    ( (c, C_elt)->

      my(idl_relnorm_hnf);
      idl_relnorm_hnf = rnfidealnormrel(shraycontextKoverK0(c), rnfidealabstorel(shraycontextKoverK0(c), idealhnf(shraycontextK(c), C_elt)));

      my(idl_relnorm_hnf_dlog, idl_relnorm_hnf_gen);
      [idl_relnorm_hnf_dlog, idl_relnorm_hnf_gen] = bnrisprincipal(shraycontextK0plus(c), idl_relnorm_hnf);

      if(idl_relnorm_hnf_dlog != 0, error("relative norm is not principal"));
      return([C_elt, idl_relnorm_hnf_gen]); /* bnrisprincipal gives positive generators so we cool*/
    )
  );
}

fgaginitshimuray(whatever, m, flag) = {

  if(m == 0, error("fill in a value for m"));

  my(cm);
  cm = cm_from_whatever(whatever);

  my(K, K0, KoverK0);
  K = cm_K(cm);
  K0 = cm_K0(cm);
  KoverK0 = rnfinit(K0, y^2 - z);
  nfinit(KoverK0);

  my(OKstar);
  OKstar = fgaginit(K, 2);
  
  my(OKmstar);
  OKmstar = fgaginit([K, m], 3);

  my(H_OKm1star);
  H_OKm1star = fgagmorker(OKstar, OKmstar, (x)->x);

  my(OK0star);
  OK0star = fgaginit(K0, 2);

  my(H_relnorm_imageof_OK1star);
  H_relnorm_imageof_OK1star = fgagmorsubgpimg(
    OKstar,
    OK0star,
    (x)->rnfeltnorm( KoverK0, x ),
    H_OKm1star
  );

  my(H_OK0plusstar);
  my(eps0, embsign);
  
  eps0 = K0.fu[1];
  embsign = embedsign(K0, eps0);
  
  H_OK0plusstar = matdiagonal([2, 2]);
  if( embsign > 0, H_OK0plusstar = matdiagonal([2, 1]) );
  if( embsign < 0, H_OK0plusstar = [2, 1; 0, 1] );

  my(coker1);
  coker1 = fgagsubgpquo(OK0star, H_relnorm_imageof_OK1star, H_OK0plusstar);

  my(Km);
  Km = bnrinit(K, m, 1);

  my(ClKm);
  ClKm = fgaginit(Km, 1);

  my(K0plus);
  K0plus = bnrinit(K0, [1, [1, 1]], 1);

  my(ClK0plus);
  ClK0plus = fgaginit(K0plus, 1);

  my(H_kerphi2);
  H_kerphi2 = fgagmorker(
    ClKm,
    ClK0plus,
    ((x)->rnfidealnormrel(KoverK0, idealhnf(K, x))) /*changed*/
  );

  my(ker2);
  ker2 = fgagsubgpsnf(H_kerphi2, ClKm);

  my(eff, gee, effinv, geeinv);
  eff = shrayeff();
  gee = shraygee();
  effinv = shrayeffinv();
  geeinv = shraygeeinv();

  my(mul, pow, con, mulpowcon);
  mul = shraymul();
  pow = shraypow();
  con = [ coker1, ker2, K, KoverK0, Km, K0plus ];
  if(flag == 1, con = [coker1, ker2, K, KoverK0, Km, K0plus, [cm, [OKstar, OKmstar, OK0star], [ ClKm, ClK0plus] ] ]; );
  mulpowcon = [mul, pow, con];

  /*breakpoint();*/

  my(ret);
  ret = fgaggrpext(coker1, ker2, eff, effinv, gee, geeinv, mulpowcon, 0);
  return(ret);

}

/* ShGrp = fgaginitshimuray([5, 134, 4169], 2, 1); */

shraycontextcm(Grp) = { return((xxgrpextinputcon(xxfgagcon(Grp)))[7][1]) };
shraycontextOKstar(Grp) = { return((xxgrpextinputcon(xxfgagcon(Grp)))[7][2][1]) };
shraycontextOKmstar(Grp) = { return((xxgrpextinputcon(xxfgagcon(Grp)))[7][2][2]) };
shraycontextOK0star(Grp) = { return((xxgrpextinputcon(xxfgagcon(Grp)))[7][2][3]) };
shraycontextClKm(Grp) = { return((xxgrpextinputcon(xxfgagcon(Grp)))[7][3][1]) };
shraycontextClK0plus(Grp) = { return((xxgrpextinputcon(xxfgagcon(Grp)))[7][3][2]) };

shraymor(c, x) = {

  my(Kr, KroverK0r, LoverK, LoverKr, sigm);
  Kr = c[7];
  KroverK0r = c[8];
  LoverK = c[11];
  LoverKr = c[12];
  sigm = c[13];

  my(a, n);

  if (LoverKr == 0, /* Galois case */
      a = idealmul (Kr, x, nfgaloisapply (Kr, sigm, x)),
   /* else non Galois case */
      a = rnfidealnormrel (LoverK, rnfidealup (LoverKr, idealhnf(Kr, x)));
   );
   n = idealnorm (Kr, x);

 return([a, n]);

}

/*
  V can be:
    a "t_VEC" with two elements : [A, B]
    a "t_VEC" with three elements : [D, A, B]
  m can be:
    a "t_INT": an integer
    a "t_MAT": the HNF of an ideal of K
    a "t_VEC": the Shimura ray class group of [D, A, B] mod m
*/
shfldcontainshcf(V, m, flag) = {

  my(ShGrp);
  if(#V == 5,
    ShGrp = V,
    ShGrp = fgaginitshimuray(V, m, 1);
  );

  my(cm, K0r, Kr, KroverK0r, Krm, K0rm);
  cm = shraycontextcm(ShGrp);
  K0r = cm[6];
  Kr = cm[7];
  KroverK0r = cm[8];
  nfinit(KroverK0r);
  Krm = iferr(bnrinit(Kr, m, 1), EE, Kr = bnfinit(Kr.pol, 1); bnrinit(Kr, m, 1));
  K0rm = bnrinit(K0r, m, 1);

  my(ClKr, ClKrm, ClK0rm);
  ClKr = fgaginit(Kr, 1);
  ClKrm = fgaginit(Krm, 1);
  ClK0rm = fgaginit(K0rm, 1);

  my(H_TypNmSubGp, H_HK0rm, H_HK1, H_int);
  H_TypNmSubGp = fgagmorker( ClKrm, ShGrp, ((x)->shraymor(cm,x)) );
  H_HK0rm = fgagmorker( ClKrm, ClK0rm, ((x)->rnfidealnormrel(KroverK0r, idealhnf(Kr, x) )) );
  H_HK1 = fgagmorker( ClKrm, ClKr, ((x)->x) );
  H_int = fgagsubgpint(ClKrm, H_TypNmSubGp, H_HK0rm);

  my(ret = fgagissubgp(ClKrm, H_int, H_HK1));
  if(flag == 1, return( [ret, V, m, ClKr, ClKrm, ClK0rm, ShGrp, H_TypNmSubGp, H_HK0rm, H_HK1, H_int] ) );
  return(ret);
}

shfldcontainshcfm0m1(V, m0, m1, flag) = {

  my(ShGrp);
  if(#V == 5,
    ShGrp = V,
    ShGrp = fgaginitshimuray(V, m1, 1);
  );

  my(cm, K0r, Kr, KroverK0r, Krm, K0rm);
  cm = shraycontextcm(ShGrp);
  K0r = cm[6];
  Kr = cm[7];
  KroverK0r = cm[8];
  nfinit(KroverK0r);
  Krm = bnrinit(Kr, max(m0,m1), 1);
  K0rm = bnrinit(K0r, m0, 1);

  my(ClKr, ClKrm, ClK0rm);
  ClKr = fgaginit(Kr, 1);
  ClKrm = fgaginit(Krm, 1);
  ClK0rm = fgaginit(K0rm, 1);

  my(H_TypNmSubGp, H_HK0rm, H_HK1, H_int);
  H_TypNmSubGp = fgagmorker( ClKrm, ShGrp, ((x)->shraymor(cm,x)) );
  H_HK0rm = fgagmorker( ClKrm, ClK0rm, ((x)->rnfidealnormrel(KroverK0r, idealhnf(Kr, x) )) );
  H_HK1 = fgagmorker( ClKrm, ClKr, ((x)->x) );
  H_int = fgagsubgpint(ClKrm, H_TypNmSubGp, H_HK0rm);

  my(ret = fgagissubgp(ClKrm, H_int, H_HK1));
  if(flag == 1, return( [ret, V, m, ClKr, ClKrm, ClK0rm, ShGrp, H_TypNmSubGp, H_HK0rm, H_HK1, H_int] ) );
  return(ret);
}

cmkcontainshcf(V, m, flag) = {

  my(ShGrp);
  if(#V == 5,
    ShGrp = V,
    ShGrp = fgaginitshimuray(V, m, 1);
  );

  my(cm, K0r, Kr, KroverK0r, Krm, K0rm);
  cm = shraycontextcm(ShGrp);
  K0r = bnfinit(cm[6].pol, 1);
  KroverK0r = rnfinit(K0r, y^2 - z);
  Kr = bnfinit(nfinit(KroverK0r).pol, 1);
  Krm = bnrinit(Kr, m, 1);
  K0rm = bnrinit(K0r, m, 1);

  my(ClKr, ClKrm, ClK0rm);
  ClKr = fgaginit(Kr, 1);
  ClKrm = fgaginit(Krm, 1);
  ClK0rm = fgaginit(K0rm, 1);

  my(H_TypNmSubGp, H_HK0rm, H_HK1, H_int);
  H_TypNmSubGp = fgagmorker( ClKrm, ShGrp, ((x)->shraymor(cm,x)) );
  H_HK0rm = fgagmorker( ClKrm, ClK0rm, ((x)->rnfidealnormrel(KroverK0r, idealhnf(Kr, x) )) );
  H_HK1 = fgagmorker( ClKrm, ClKr, ((x)->x) );
  H_int = fgagsubgpint(ClKrm, H_TypNmSubGp, H_HK0rm);

  my(ret = fgagissubgp(ClKrm, H_TypNmSubGp, H_HK1));
  if(flag == 1, return( [ret, V, m, ClKr, ClKrm, ClK0rm, ShGrp, H_TypNmSubGp, H_HK0rm, H_HK1, H_int] ) );
  return(ret);
}


shraysubgptypenorm(ShGrp, m, flag) = {

  my(cm, K0r, Kr, KroverK0r, Krm, K0rm);
  cm = shraycontextcm(ShGrp);
  K0r = cm[6];
  Kr = cm[7];
  KroverK0r = cm[8];
  nfinit(KroverK0r);
  Krm = bnrinit(Kr, m, 1);
  K0rm = bnrinit(K0r, m, 1);

  my(ClKr, ClKrm, ClK0rm);
  ClKr = fgaginit(Kr, 1);
  ClKrm = fgaginit(Krm, 1);
  ClK0rm = fgaginit(K0rm, 1);

  my(H_TypNmSubGp, H_HK0rm, H_HK1, H_int);
  H_TypNmSubGp = fgagmorker( ClKrm, ShGrp, ((x)->shraymor(cm,x)) );
  if(flag == 0, return([H_TypNmSubGp, ClKrm]) );
  return( fgagsubgpsnf(H_TypNmSubGp, ClKrm) );

}

install(idealmoddivisor, GG);

shraycosetsmacronaive(ShGrp, m) = {

  my(cm, K0r, Kr, KroverK0r, Krm, K0rm);
  cm = shraycontextcm(ShGrp);
  K0r = cm[6];
  Kr = cm[7];
  KroverK0r = cm[8];
  nfinit(KroverK0r);
  Krm = bnrinit(Kr, m, 1);
  K0rm = bnrinit(K0r, m, 1);

  my(ClKr, ClKrm, ClK0rm);
  ClKr = fgaginit(Kr, 1);
  ClKrm = fgaginit(Krm, 1);
  ClK0rm = fgaginit(K0rm, 1);

  my(H_TypNmSubGp, H_HK0rm, H_HK1, H_int);
  H_TypNmSubGp = fgagmorker( ClKrm, ShGrp, ((x)->shraymor(cm,x)) );
  H_HK0rm = fgagmorker( ClKrm, ClK0rm, ((x)->rnfidealnormrel(KroverK0r, idealhnf(Kr, x) )) );
  H_HK1 = fgagmorker( ClKrm, ClKr, ((x)->x) );
  H_int = fgagsubgpint(ClKrm, H_TypNmSubGp, H_HK0rm);

  my(ClKrmquo);
  ClKrmquo = fgagsubgpquo(ClKrm, H_TypNmSubGp, matid(#H_TypNmSubGp));

  my(H_HK1_subgpofkercalN);
  H_HK1_subgpofkercalN = fgagmorker( ClKrmquo, ClKr, ((x)->x) );

  my(unreduced);
  unreduced = fgaglistcoseltsnaive(ClKrmquo, ClKrm, H_HK1_subgpofkercalN, 1);

  return( apply(xx->shraycosetreduce(ClKrm, xx), unreduced) );
}

shraycosetreduce(ClKrm, elts) = {

  my(dlooog, idealstoreduce);
  [dlooog, idealstoreduce] = elts;

  my(Krm);
  Krm = ClKrm[5];

  my(reducedideals);
  vector(#idealstoreduce, i, [idealstoreduce[i][1], idealmoddivisor(Krm, idealstoreduce[i][2])] );

  return( [dlooog, idealstoreduce] );
}

/* FAILS ON HilbertClassFieldusingCM([28, 12, 29], bpind=7, intermediate=True); */
isprobablygalois(K, pol, plist, p_ub) = {

  if(p_ub == 0, p_ub = 100);

  one = 0;
  forprime(p=2, p_ub,
    if( K.disc%p==0, next() );
    if( setsearch(plist, p)!=0, next() );
    id=idealprimedec( K, p );
    for(j=1,#id,
      my(F); /* F exists as global variable in SAGE script! */
      F = nffactormod(K, pol, id[j]);
      tmp = vecsort(apply(poldegree,F[,1]),,8);
      if(#tmp != 1, return(Strprintf("POLYNOMIAL FOUND DOES NOT SEEM TO BE GALOIS! Prime %s", p)));
      if(vecsort(F[,2])[1] == 1, one = 1);
    );
  );
  if(one != 1, return(Strprintf("POLYNOMIAL DOES NOT SEEM TO BE IRREDUCIBLE!")));
  return(1);

}