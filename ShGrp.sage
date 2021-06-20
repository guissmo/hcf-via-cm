gp.read(FGAGDIR+"/fgag.gp");
gp.read(FGAGDIR+"/fgagshimuray.gp");

class ShGrp:

    def __init__(self, V, m):
        self.V = V;
        self.m = m;
        self.gp = gp.fgaginitshimuray(V, m, 1);
        [self.coker1, self.coker2, self.K, self.KoverK0, self.Km, self.K0plus, [self.cm, [self.OKstar, self.OKmstar, self.OK0star], [self.ClKm, self.ClK0plus] ]] = self.gp[5][7]

    def __repr__(self):
        return(str(f"Shimura ray class group of CM field defined by {str(self.K[7][1])} for the modulus {str(self.m)}"))

class ShGrps:
    """
    |
    |
    """

    def __init__(self, V):
        self.V = V
        self.cm = gp.cm_from_whatever(V)
        self.K0r1 = self.cm[6]
        self.Kr1 = self.cm[7]

        self.ShGrpm = {}
        self.K0rm = {}
        self.Krm = {}
        self.ClK0rm = {}
        self.ClKrm = {}

    def ShGrp(self, m=1):
        if m not in self.ShGrpm:
            self.ShGrpm[m] = ShGrp(self.cm, m)
        return(self.ShGrpm[m])

    def K0r(self, m=1):
        if m not in self.K0rm:
            self.K0rm[m] = gp.bnrinit(self.K0r1, m, 1)
        return(self.K0rm[m])

    def Kr(self, m=1):
        if m not in self.Krm:
            self.Krm[m] = gp.bnrinit(self.Kr1, m, 1)
        return(self.Krm[m])

    def ClK0r(self, m=1):
        if m not in self.ClK0rm:
            self.ClK0rm[m] = gp.fgaginit(self.K0r(m), 1)
        return(self.ClK0rm[m])

    def ClKr(self, m=1):
        if m not in self.ClKrm:
            self.ClKrm[m] = gp.fgaginit(self.Kr(m), 1)
        return(self.ClKrm[m])

    def typenormsubgroup_hnf(self, m0, m1):
        # print(f"((x)->shraymor({str(self.ShGrp(m1).cm)},x))")
        return( gp.fgagmorker(self.ClKr(m0), self.ShGrp(m1).gp, gp(f"((x)->shraymor({str(gp.substpol(self.ShGrp(m1).cm,gp('y'),gp('z')))},x))")) )

        
