from fractions import Fraction
import copy

class Raw:
    def __init__(self, num: int, den: int, real: bool):
        if den != 0:
            self.num = num
            self.den = den
            self.real = real
        else:
            raise ZeroDivisionError
        self.simplify()

    def simplify(self):
        temp1 = Fraction(self.num, self.den).numerator
        temp2 = Fraction(self.num, self.den).denominator
        self.num, self.den = temp1, temp2

    def __str__(self):
        if self.den == 1 and self.real:
            return str(self.num)
        elif self.den == 1 and not self.real:
            return str(self.num) + "i"
        elif self.den != 1 and self.real:
            return str(self.num) + "/" + str(self.den)
        elif self.den != 1 and not self.real:
            return str(self.num) + "i" + "/" + str(self.den)
        else:
            raise Exception("UnknownError")
    __repr__ = __str__

    def __mul__(self, other):
        fra = Fraction(self.num, self.den) * Fraction(other.num, other.den)
        ree = self.real + other.real
        if ree == 2:
            return Raw(fra.numerator, fra.denominator, True)
        elif ree == 1:
            return Raw(fra.numerator, fra.denominator, False)
        elif ree == 0:
            return Raw(-fra.numerator, fra.denominator, True)
        else:
            raise Exception("UnknownError")

    def __eq__(self, other):
        temp = (self.num == other.num and self.den == other.den and self.real == other.real)
        return (temp or self.num == other.num == 0)

class Number:
    def __init__(self, *args: Raw, simp=False):
        self.args = [*args]
        self.simp = simp
        if not simp:
            self.args = self.simplify().args

    def simplify(self):
        tempT = sum([Fraction(i.num, i.den) for i in self.args if i.real])
        zeroT = ((tempT == 0) or ((type(tempT) == Fraction) and (tempT.numerator == 0)))
        if not zeroT:
            tempT = Raw(tempT.numerator, tempT.denominator, True)

        tempF = sum([Fraction(i.num, i.den) for i in self.args if not i.real])
        zeroF = ((tempF == 0) or ((type(tempF) == Fraction) and (tempF.numerator == 0)))
        if not zeroF:
            tempF = Raw(tempF.numerator, tempF.denominator, False)

        if zeroT and zeroF:
            return Number(Raw(0,1,True), simp=True)
        elif zeroF:
            return Number(tempT, simp=True)
        elif zeroT:
            return Number(tempF, simp=True)
        else:
            return Number(tempT, tempF, simp=True)

    def __str__(self):
        if len(self.args) == 0:
            return ""
        elif len(self.args) == 1:
            return self.args[0].__str__()
        elif len(self.args) == 2:
            return "(" + self.args[0].__str__() + "+" + self.args[1].__str__() + ")"
        else:
            raise Exception("UnknownError")
    __repr__ = __str__

    def __mul__(self, other):
        temp1 = self.args
        temp2 = other.args
        return Number(*[i * j for i in temp1 for j in temp2])

    def smul(self, other):
        if type(other) == Raw:
            temp = Number(*self.args, simp=True)
            for i,j in enumerate(temp.args):
                temp.args[i] = temp.args[i] * other
            return temp
        elif type(other) == int:
            return self.smul(Raw(other,1,True))
        else:
            raise Exception("UnknownError")

    def __eq__(self, other):
        return self.args == other.args

    def __add__(self, other):
        return Number(*(self.args + other.args))

    def __copy__(self):
        temp = []
        for i in self.args:
            temp = temp + [Raw(i.num, i.den, i.real)]
        return Number(*temp)

class Symbol:
    def __init__(self, name: str, coef = Number(Raw(1,1,True)),
                 ind_s = 0, ind_sp = 0, ind_g = (), ind_gp = (),
                 dis = "?"):
        self.name = name
        self.coef = coef
        self.ind_s = ind_s
        self.ind_sp = ind_sp
        self.ind_g = ind_g
        self.ind_gp = ind_gp
        self.dis = dis

    def smul(self, other):
        if type(other) == Number:
            temp = Symbol(self.name, coef=self.coef,
                          ind_s=self.ind_s, ind_sp=self.ind_sp,
                          ind_g=self.ind_g, ind_gp=self.ind_gp, dis=self.dis)
            temp.coef = temp.coef * other
            return temp
        elif type(other) == Raw:
            return self.smul(Number(other))
        elif type(other) == int:
            return self.smul(Raw(other,1,True))
        else:
            raise Exception("UnknownError")

    def __str__(self):
        return str(self.coef) + " " + self.dis
    __repr__ = __str__

    def __copy__(self):
        return Symbol(self.name, copy.copy(self.coef),
                      self.ind_s, self.ind_sp, self.ind_g, self.ind_gp, self.dis)

    def __eq__(self, other):
        return (self.name == other.name and self.coef == other.coef
                and self.ind_s == other.ind_s and self.ind_sp == other.ind_sp
                and self.ind_g == other.ind_g and self.ind_gp == other.ind_gp
                and self.dis == other.dis)

one_sym = Symbol("one", ind_s = 0, ind_g = (0,), dis = "")

delta = Symbol("delta", ind_s = 0, ind_g = (2,), dis = "\u03B4")
tiger = Symbol("tiger", ind_s = 0, ind_g = (3,), dis = "f")

met = Symbol("met", ind_s = 2, ind_g = (0,), dis = "\u03B7")
levi = Symbol("levi", ind_s = 4, ind_g = (0,), dis = "\u03B5")

catn = Symbol("catn", ind_s = 4, ind_g = (0,), dis = "cat[\u03B7\u03B7]")
catp = Symbol("catp", ind_s = 4, ind_g = (0,), dis = "(cat+)")
catm = Symbol("catm", ind_s = 4, ind_g = (0,), dis = "(cat-)")

dogn = Symbol("dogn", ind_s = 4, ind_g = (0,), dis = "(dog)")
dogp = Symbol("dogp", ind_s = 4, ind_g = (0,), dis = "(dog+)")
dogm = Symbol("dogm", ind_s = 4, ind_g = (0,), dis = "(dog-)")

foxn = Symbol("foxn", ind_s = 4, ind_g = (0,), dis = "(fox)")
foxp = Symbol("foxp", ind_s = 4, ind_g = (0,), dis = "(fox+)")
foxm = Symbol("foxm", ind_s = 4, ind_g = (0,), dis = "(fox-)")

eel = Symbol("eel", ind_s = 6, ind_g = (0,), dis = "eel")
elk = Symbol("elk", ind_s = 6, ind_g = (0,), dis = "elk")
owlp = Symbol("owlp", ind_s = 6, ind_g = (0,), dis = "(owl+)")
owlm = Symbol("owlm", ind_s = 6, ind_g = (0,), dis = "(owl-)")

class Gamma:
    def __init__(self,
                 c1 = Number(Raw(0,1,True)), c2 = Number(Raw(0,1,True)),
                 c3 = Number(Raw(0,1,True)), c4 = Number(Raw(0,1,True)),
                 c5 = Number(Raw(0,1,True)), c6 = Number(Raw(0,1,True))):
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3
        self.c4 = c4
        self.c5 = c5
        self.c6 = c6

    def __str__(self):
        temp = ""
        if self.c1 != Number(Raw(0,1,True)):
            if temp == "" and self.c1 != Number(Raw(1,1,True)):
                temp += f"{self.c1} \u03B3sl"
            elif self.c1 != Number(Raw(1,1,True)):
                temp += f" + {self.c1} \u03B3sl"
            elif temp == "":
                temp += f"\u03B3sl"
            else:
                temp += f" + \u03B3sl"
        if self.c2 != Number(Raw(0,1,True)):
            if temp == "" and self.c2 != Number(Raw(1,1,True)):
                temp += f"{self.c2} \u03B3sr"
            elif self.c2 != Number(Raw(1,1,True)):
                temp += f" + {self.c2} \u03B3sr"
            elif temp == "":
                temp += f"\u03B3sr"
            else:
                temp += f" + \u03B3sr"
        if self.c3 != Number(Raw(0,1,True)):
            if temp == "" and self.c3 != Number(Raw(1,1,True)):
                temp += f"{self.c3} \u03B3vl"
            elif self.c3 != Number(Raw(1,1,True)):
                temp += f" + {self.c3} \u03B3vl"
            elif temp == "":
                temp += f"\u03B3vl"
            else:
                temp += f" + \u03B3vl"
        if self.c4 != Number(Raw(0,1,True)):
            if temp == "" and self.c4 != Number(Raw(1,1,True)):
                temp += f"{self.c4} \u03B3vr"
            elif self.c4 != Number(Raw(1,1,True)):
                temp += f" + {self.c4} \u03B3vr"
            elif temp == "":
                temp += f"\u03B3vr"
            else:
                temp += f" + \u03B3vr"
        if self.c5 != Number(Raw(0,1,True)):
            if temp == "" and self.c5 != Number(Raw(1,1,True)):
                temp += f"{self.c5} \u03B3tl"
            elif self.c5 != Number(Raw(1,1,True)):
                temp += f" + {self.c5} \u03B3tl"
            elif temp == "":
                temp += f"\u03B3tl"
            else:
                temp += f" + \u03B3tl"
        if self.c6 != Number(Raw(0,1,True)):
            if temp == "" and self.c6 != Number(Raw(1,1,True)):
                temp += f"{self.c6} \u03B3tr"
            elif self.c6 != Number(Raw(1,1,True)):
                temp += f" + {self.c6} \u03B3tr"
            elif temp == "":
                temp += f"\u03B3tr"
            else:
                temp += f" + \u03B3tr"
        return temp
    __repr__ = __str__

    def __add__(self, other):
        t1 = Number(*(self.c1.args + other.c1.args))
        t2 = Number(*(self.c2.args + other.c2.args))
        t3 = Number(*(self.c3.args + other.c3.args))
        t4 = Number(*(self.c4.args + other.c4.args))
        t5 = Number(*(self.c5.args + other.c5.args))
        t6 = Number(*(self.c6.args + other.c6.args))
        return Gamma(t1, t2, t3, t4, t5, t6)

    def trimul(self, self2, self3):
        t_one = (self.c1 * self2.c1 * self3.c1).smul(2)\
             + (self.c2 * self2.c2 * self3.c2).smul(2)
        T_one = Term(t_one, [], [], [])

        t_met = (self.c2 * self2.c3 * self3.c4).smul(2)\
             + (self.c4 * self2.c2 * self3.c3).smul(2)\
             + (self.c3 * self2.c4 * self3.c2).smul(2)\
             + (self.c1 * self2.c4 * self3.c3).smul(2)\
             + (self.c3 * self2.c1 * self3.c4).smul(2)\
             + (self.c4 * self2.c3 * self3.c1).smul(2)
        T_met = Term(t_met, [], [met], [])

        t_catp = (self.c1 * self2.c5 * self3.c5).smul(2)\
                 + (self.c5 * self2.c1 * self3.c5).smul(2)\
                 + (self.c5 * self2.c5 * self3.c1).smul(2)\
                 + (self.c5 * self2.c4 * self3.c3).smul(Raw(-2,1,False))\
                 + (self.c4 * self2.c3 * self3.c5).smul(Raw(-2,1,False))
        T_catp = Term(t_catp, [], [catp], [])

        t_catm = (self.c2 * self2.c6 * self3.c6).smul(2)\
                 + (self.c6 * self2.c2 * self3.c6).smul(2)\
                 + (self.c6 * self2.c6 * self3.c2).smul(2)\
                 + (self.c6 * self2.c3 * self3.c4).smul(Raw(-2,1,False))\
                 + (self.c3 * self2.c4 * self3.c5).smul(Raw(-2,1,False))
        T_catm = Term(t_catm, [], [catm], [])

        t_dogp = (self.c4 * self2.c6 * self3.c3).smul(Raw(-2,1,False))
        T_dogp = Term(t_dogp, [], [dogp], [])

        t_dogm = (self.c3 * self2.c5 * self3.c4).smul(Raw(-2,1,False))
        T_dogm = Term(t_dogm, [], [dogm], [])

        t_owlp = (self.c5 * self2.c5 * self3.c5).smul(2)
        T_owlp = Term(t_owlp, [], [owlp], [])

        t_owlm = (self.c6 * self2.c6 * self3.c6).smul(2)
        T_owlm = Term(t_owlm, [], [owlm], [])

        temp = [[t_one, copy.copy(one_sym)], [t_met, copy.copy(met)],
                [t_catp, copy.copy(catp)], [t_catm, copy.copy(catm)],
                [t_dogp, copy.copy(dogp)], [t_dogm, copy.copy(dogm)],
                [t_owlp, copy.copy(owlp)], [t_owlm, copy.copy(owlm)]]
        return [i for i in temp if not (i[0] == Number(Raw(0, 1, True)))]
        #return Expression(T_one, T_met, T_catp, T_catm, T_dogp, T_dogm, T_owlp, T_owlm)

class Gener:
    def __init__(self,
                 c1 = Number(Raw(0,1,True)), c2 = Number(Raw(0,1,True))):
        self.c1 = c1
        self.c2 = c2

    def __str__(self):
        temp = ""
        if self.c1 != Number(Raw(0,1,True)):
            if temp == "" and self.c1 != Number(Raw(1,1,True)):
                temp += f"{self.c1} \u03B4"
            elif self.c1 != Number(Raw(1,1,True)):
                temp += f" + {self.c1} \u03B4"
            elif temp == "":
                temp += f"\u03B4"
            else:
                temp += f" + \u03B4"
        if self.c2 != Number(Raw(0,1,True)):
            if temp == "" and self.c2 != Number(Raw(1,1,True)):
                temp += f"{self.c2} \u03C4"
            elif self.c2 != Number(Raw(1,1,True)):
                temp += f" + {self.c2} \u03C4"
            elif temp == "":
                temp += f"\u03C4"
            else:
                temp += f" + \u03C4"
        return temp
    __repr__ = __str__

    def __add__(self, other):
        t1 = Number(*(self.c1.args + other.c1.args))
        t2 = Number(*(self.c2.args + other.c2.args))
        return Gamma(t1, t2)

    def trimul(self, self2, self3):
        t_one = (self.c1 * self2.c1 * self3.c1).smul(2)

        t_delta = (self.c1 * self2.c2 * self3.c2).smul(Raw(1,2,True)) \
                  + (self.c2 * self2.c1 * self3.c2).smul(Raw(1,2,True)) \
                  + (self.c2 * self2.c2 * self3.c1).smul(Raw(1,2,True))

        t_tiger = (self.c2 * self2.c2 * self3.c2).smul(Raw(1,4,False))

        temp = [[t_one, copy.copy(one_sym)], [t_delta, copy.copy(delta)], [t_tiger, copy.copy(tiger)]]
        return [i for i in temp if not (i[0] == Number(Raw(0, 1, True)))]

class Field:
    def __init__(self, name: str, kind = "c", der = 0,
                 coef = Number(Raw(1,1,True)),
                 ind_s = 0, ind_sp = 0, ind_g = (), ind_gp = ()):
        self.name = name
        self.kind = kind
        self.der = der
        self.coef = coef
        self.ind_s = ind_s
        self.ind_sp = ind_sp
        self.ind_g = ind_g
        self.ind_gp = ind_gp

    def smul(self, other):
        if type(other) == Number:
            temp = Field(self.name, kind=self.kind, der=self.der,
                         coef=self.coef,
                         ind_s=self.ind_s, ind_sp=self.ind_sp,
                         ind_g=self.ind_g, ind_gp=self.ind_gp)
            temp.coef = temp.coef * other
            return temp
        elif type(other) == Raw:
            return self.smul(Number(other))
        elif type(other) == int:
            return self.smul(Raw(other,1,True))
        else:
            raise Exception("UnknownError")

    def __str__(self):
        if self.der == 0:
            return str(self.coef) + " " + self.name
        elif self.der == 1:
            return str(self.coef) + " \u2202" + self.name
        else:
            raise Exception("UnknownError")
    __repr__ = __str__

    def __lt__(self, other):
        """Canonical ordering"""
        if self.kind == other.kind == "c":
            if self.der > other.der:
                return True
            elif self.der < other.der:
                return False
            elif self.ind_s < other.ind_s:
                return True
            elif self.ind_s > other.ind_s:
                return False
            elif self.ind_sp < other.ind_sp:
                return True
            elif self.ind_sp > other.ind_sp:
                return False
            elif self.ind_g < other.ind_g:
                return True
            elif self.ind_g > other.ind_g:
                return False
            elif self.ind_gp < other.ind_gp:
                return True
            elif self.ind_gp > other.ind_gp:
                return False
            elif self.name < other.name:
                return True
            elif self.name > other.name:
                return False
            else:
                return False
        else:
            return False

    def __copy__(self):
        return Field(self.name, self.kind, self.der,
                     copy.copy(self.coef),
                     self.ind_s, self.ind_sp, self.ind_g, self.ind_gp)

    def __eq__(self, other):
        temp = (self.name == other.name and self.kind == other.kind
                and self.der == other.der and self.coef == other.coef
                and self.ind_s == other.ind_s and self.ind_sp == other.ind_sp
                and self.ind_g == other.ind_g and self.ind_gp == other.ind_gp)
        return temp

class Constant:
    pass

class Term:
    def __init__(self, nmbr: Number, consts: list,
                 syms: list, fields: list):
        self.nmbr = nmbr
        self.consts = consts
        self.syms = syms
        self.fields = fields

    def __str__(self):
        coef = copy.copy(self.nmbr)
        for i in self.consts:
            coef = coef * i.coef
        for i in self.syms:
            coef = coef * i.coef
        for i in self.fields:
            coef = coef * i.coef

        if coef != Number(Raw(1,1,True)):
            temp1 = str(coef)
        else:
            temp1 = ""
        temp2 = ""
        for i,j in enumerate(self.consts):
            if i == 0:
                temp2 = temp2 + j.name
            else:
                temp2 = temp2 + " " + j.name
        if temp1 != "" and temp2 != "":
            temp2 = " " + temp2
        temp3 = ""
        for i,j in enumerate(self.syms):
            if i == 0:
                temp3 = temp3 + j.dis
            else:
                temp3 = temp3 + " " + j.dis
        if (temp1 + temp2) != "" and temp3 != "":
            temp3 = " " + temp3
        temp4 = ""
        for i,j in enumerate(self.fields):
            if i == 0:
                temp4 = temp4 + "\u2202"*j.der + j.name
            else:
                temp4 = temp4 + " " + "\u2202"*j.der + j.name
        if (temp1 + temp2 + temp3) != "" and temp4 != "":
            temp4 = " " + temp4
        return temp1 + temp2 + temp3 + temp4
    __repr__ = __str__

    def __copy__(self):
        return Term(copy.copy(self.nmbr), copy.deepcopy(self.consts),
                    copy.deepcopy(self.syms), copy.deepcopy(self.fields))

    def __eq__(self, other):
        """Linear dependence"""
        return (self.consts == other.consts
                and self.syms == other.syms and self.fields == other.fields)

class Expression:
    def __init__(self, *terms: Term, simp=0):
        self.terms = [*terms]
        self.simp = simp
        if not simp == 3:
            self.terms = self.simplify().terms

    def __str__(self):
        temp = ""
        for i,j in enumerate(self.terms):
            if i == 0:
                temp = temp + str(j)
            else:
                temp = temp + " + " + str(j)
        return str(temp)
    __repr__ = __str__

    def simplify(self):
        simp0 = self.simp
        terms = self.terms.copy()

        """Collect the coefficients to the number in front"""
        for i in terms:
            for j in i.consts:
                if not j.coef == Number(Raw(1, 1, True)):
                    i.nmbr = i.nmbr * j.coef
                    j.coef = Number(Raw(1, 1, True))
            for j in i.syms:
                if not j.coef == Number(Raw(1, 1, True)):
                    i.nmbr = i.nmbr * j.coef
                    j.coef = Number(Raw(1, 1, True))
            for j in i.fields:
                if not j.coef == Number(Raw(1, 1, True)):
                    i.nmbr = i.nmbr * j.coef
                    j.coef = Number(Raw(1, 1, True))

        """Remove terms with number 0"""
        temp = []
        for i in terms:
            if i.nmbr != Number(Raw(0,1,True)):
                temp.append(i)
        terms = temp

        """Remove terms with derivatives only"""
        temp = []
        for i in terms:
            if len(i.fields) == 0:
                temp.append(i)
            else:
                t = 0
                for j in i.fields:
                    if j.kind == "d":
                        t += 1
                if t < len(i.fields):
                    temp.append(i)
        terms = temp

        """Reorder the terms with two derivatives, one field"""
        temp = []
        for i in terms:
            if len(i.fields) == 3 and len([j for j in i.fields if j.kind == "d"]) == 2:
                if i.fields[0].kind != "d":
                    i.nmbr = i.nmbr.smul(Raw(1, 2, True))
                    t1 = copy.copy(i)
                    tempa = t1.fields[0].ind_s
                    t1.fields = [t1.fields[1], t1.fields[2], t1.fields[0]]
                    t1.syms[0] = PermRule(t1.syms[0], Permutation(tempa + 2, tuple(
                        [2] + list(range(3, tempa + 3)) + [1]
                    ))).evaluate()
                    temp.append(t1)

                    t2 = copy.copy(i)
                    tempa = t2.fields[0].ind_s
                    t2.fields = [t2.fields[2], t2.fields[1], t2.fields[0]]
                    t2.syms[0] = PermRule(t2.syms[0], Permutation(tempa + 2, tuple(
                        list(range(3, tempa + 3)) + [2] + [1]
                    ))).evaluate()
                    temp.append(t2)

                elif i.fields[1].kind != "d":
                    i.nmbr = i.nmbr.smul(Raw(1, 2, True))
                    t1 = copy.copy(i)
                    tempa = t1.fields[1].ind_s
                    t1.fields = [t1.fields[2], t1.fields[0], t1.fields[1]]
                    t1.syms[0] = PermRule(t1.syms[0], Permutation(tempa + 2, tuple(
                        list(range(3, tempa + 3)) + [1] + [2]
                    ))).evaluate()
                    temp.append(t1)

                    t2 = copy.copy(i)
                    tempa = t2.fields[1].ind_s
                    t2.fields = [t2.fields[0], t2.fields[2], t2.fields[1]]
                    t2.syms[0] = PermRule(t2.syms[0], Permutation(tempa + 2, tuple(
                        [1] + list(range(3, tempa + 3)) + [2]
                    ))).evaluate()
                    temp.append(t2)

                elif i.fields[2].kind != "d":
                    i.nmbr = i.nmbr.smul(Raw(1, 2, True))
                    t1 = copy.copy(i)
                    tempa = t1.fields[2].ind_s
                    t1.fields = [t1.fields[0], t1.fields[1], t1.fields[2]]
                    t1.syms[0] = PermRule(t1.syms[0], Permutation(tempa + 2, tuple(
                        [1] + [2] + list(range(3, tempa + 3))
                    ))).evaluate()
                    temp.append(t1)

                    t2 = copy.copy(i)
                    tempa = t2.fields[2].ind_s
                    t2.fields = [t2.fields[1], t2.fields[0], t2.fields[2]]
                    t2.syms[0] = PermRule(t2.syms[0], Permutation(tempa + 2, tuple(
                        [2] + [1] + list(range(3, tempa + 3))
                    ))).evaluate()
                    temp.append(t2)
                else:
                    raise Exception("UnknownError")
            else:
                temp.append(i)
        terms = temp

        """What to do about the terms with one derivative, one field"""
        for i in terms:
            if len(i.fields) == 2 and len([j for j in i.fields if j.kind == "d"]) == 1:
                raise Exception("1D1F-NotImplemented")

        """Reorder the terms with one derivative, two fields"""
        temp = []
        for i in terms:
            if len(i.fields) == 3 and len([j for j in i.fields if j.kind == "d"]) == 1:
                if i.fields[0].kind == "d":
                    if i.fields[1] < i.fields[2]:
                        tempa = i.fields[1].ind_s
                        tempb = i.fields[2].ind_s
                        temga = i.fields[1].ind_g[0]
                        temgb = i.fields[2].ind_g[0]
                        i.fields = [i.fields[2], i.fields[0], i.fields[1]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + 1, tuple(
                            list(range(tempa + 2, tempa + tempb + 2))
                            + [1] + list(range(2, tempa + 2))
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb, tuple(
                            list(range(temga + 1, temga + temgb + 1))
                            + list(range(1, temga + 1))
                        ))).evaluate()
                        temp.append(i)
                    elif i.fields[2] < i.fields[1]:
                        tempa = i.fields[2].ind_s
                        tempb = i.fields[1].ind_s
                        temga = i.fields[2].ind_g[0]
                        temgb = i.fields[1].ind_g[0]
                        i.fields = [i.fields[1], i.fields[2], i.fields[0]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + 1, tuple(
                            list(range(2, tempb + 2))
                            + list(range(tempb + 2, tempa + tempb + 2)) + [1]
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb, tuple(
                            list(range(1, temgb + 1))
                            + list(range(temgb + 1, temga + temgb + 1))
                        ))).evaluate()
                        temp.append(i)
                    else:
                        tempa = i.fields[1].ind_s
                        tempb = i.fields[2].ind_s
                        temga = i.fields[1].ind_g[0]
                        temgb = i.fields[2].ind_g[0]
                        i.fields = [i.fields[2], i.fields[0], i.fields[1]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + 1, tuple(
                            list(range(tempa + 2, tempa + tempb + 2))
                            + [1] + list(range(2, tempa + 2))
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb, tuple(
                            list(range(temga + 1, temga + temgb + 1))
                            + list(range(1, temga + 1))
                        ))).evaluate()
                        temp.append(i)
                elif i.fields[1].kind == "d":
                    if i.fields[2] < i.fields[0]:
                        temp.append(i)
                    elif i.fields[0] < i.fields[2]:
                        tempa = i.fields[0].ind_s
                        tempb = i.fields[2].ind_s
                        temga = i.fields[0].ind_g[0]
                        temgb = i.fields[2].ind_g[0]
                        i.fields = [i.fields[2], i.fields[0], i.fields[1]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + 1, tuple(
                            list(range(tempa + 2, tempa + tempb + 2))
                            + list(range(1, tempa + 1)) + [tempa + 1]
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb, tuple(
                            list(range(temga + 1, temga + temgb + 1))
                            + list(range(1, temga + 1))
                        ))).evaluate()
                        temp.append(i)
                    else:
                        temp.append(i)
                elif i.fields[2].kind == "d":
                    if i.fields[1] < i.fields[0]:
                        temp.append(i)
                    elif i.fields[0] < i.fields[1]:
                        tempa = i.fields[0].ind_s
                        tempb = i.fields[1].ind_s
                        temga = i.fields[0].ind_g[0]
                        temgb = i.fields[1].ind_g[0]
                        i.fields = [i.fields[1], i.fields[2], i.fields[0]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + 1, tuple(
                            list(range(tempa + 1, tempa + tempb + 1))
                            + [tempa + tempb + 1] + list(range(1, tempa + 1))
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb, tuple(
                            list(range(temga + 1, temga + temgb + 1))
                            + list(range(1, temga + 1))
                        ))).evaluate()
                        temp.append(i)
                    else:
                        tempa = i.fields[0].ind_s
                        tempb = i.fields[1].ind_s
                        temga = i.fields[0].ind_g[0]
                        temgb = i.fields[1].ind_g[0]
                        i.fields = [i.fields[1], i.fields[2], i.fields[0]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + 1, tuple(
                            list(range(tempa + 1, tempa + tempb + 1))
                            + [tempa + tempb + 1] + list(range(1, tempa + 1))
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb, tuple(
                            list(range(temga + 1, temga + temgb + 1))
                            + list(range(1, temga + 1))
                        ))).evaluate()
                        temp.append(i)
                else:
                    raise Exception("UnknownError")
            else:
                temp.append(i)
        terms = temp

        """Absorb the derivatives into the fields"""
        temp = []
        for i in terms:
            f_num = len(i.fields)
            d_pos = tuple([j for j,k in enumerate(i.fields) if k.kind == "d"])
            d_num = len(d_pos)
            if d_num == 0:
                temp.append(i)
            elif not(f_num == 3 and d_num == 1 and d_pos in [(1,), (2,)]):
                temp.append(i)
            else:
                if d_pos == (2,):
                    temp.append(i)
                else:
                    td = copy.copy(i)
                    td.fields = [td.fields[0], td.fields[2]]
                    td.fields[1].der += 1
                    td.fields[1].ind_s += 1
                    temp.append(td)

                    tr = copy.copy(i)
                    tempb = tr.fields[0].ind_s
                    tempa = tr.fields[2].ind_s
                    tr.fields = [tr.fields[0], tr.fields[2], tr.fields[1]]
                    tr.syms[0] = PermRule(tr.syms[0], Permutation(tempa + tempb + 1, tuple(
                        list(range(1, tempb + 1))
                        + list(range(tempb + 2, tempb + tempa + 2)) + [tempb + 1]
                    ))).evaluate()
                    temp.append(tr)
        terms = temp

        """Skew-symmetrize over the 2-form indices"""
        temp = []
        for i in terms:
            temp_list = [copy.copy(i)]
            f_num = len(i.fields)
            for j in range(f_num):
                if i.fields[j].ind_s == 2 and i.fields[j].der == 0:
                    tempc = sum([i.fields[ii].ind_s for ii in range(j)])
                    tempa = sum([i.fields[ii].ind_s for ii in range(j+1, f_num)])
                    temp_list_2 = []
                    for k in temp_list:
                        k.nmbr = k.nmbr.smul(Raw(1,2,True))
                        temp_list_2.append(copy.copy(k))
                        k.nmbr = k.nmbr.smul(-1)
                        k.syms[0] = PermRule(k.syms[0], Permutation(tempc + 2 + tempa, tuple(
                            list(range(1, tempc + 1)) + [tempc + 2, tempc + 1]
                            + list(range(tempc + 3, tempc + 3 + tempa))
                        ))).evaluate()
                        temp_list_2.append(copy.copy(k))
                    temp_list = temp_list_2
            temp = temp + temp_list
        terms = temp

        """Reorder the terms without derivatives"""
        temp = []
        for i in terms:
            f_num = len(i.fields)
            d_num = len([j for j in i.fields if j.kind == "d"])
            if d_num != 0 or f_num <= 1:
                temp.append(i)
            elif f_num == 2:
                if not i.fields[0] < i.fields[1]:
                    temp.append(i)
                else:
                    tempa = i.fields[0].ind_s
                    tempb = i.fields[1].ind_s
                    temga = i.fields[0].ind_g[0]
                    temgb = i.fields[1].ind_g[0]
                    i.fields = [i.fields[1], i.fields[0]]
                    i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb, tuple(
                        list(range(tempa + 1, tempa + tempb + 1)) + list(range(1, tempa + 1))
                    ))).evaluate()
                    i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb, tuple(
                        list(range(temga + 1, temga + temgb + 1)) + list(range(1, temga + 1))
                    ))).evaluate()
                    temp.append(i)
            elif f_num == 3:
                if not(i.fields[0].kind == i.fields[1].kind == i.fields[2].kind == "c"):
                    raise Exception("ImplementationError")
                not_ordered = (i.fields[0] < i.fields[1] or i.fields[1] < i.fields[2]
                               or i.fields[0] < i.fields[2])
                iter = 0
                while not_ordered and iter < 4:
                    if i.fields[0] < i.fields[1]:
                        tempa = i.fields[0].ind_s
                        tempb = i.fields[1].ind_s
                        tempx = i.fields[2].ind_s
                        temga = i.fields[0].ind_g[0]
                        temgb = i.fields[1].ind_g[0]
                        temgx = i.fields[2].ind_g[0]
                        i.fields = [i.fields[1], i.fields[0], i.fields[2]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + tempx, tuple(
                            list(range(tempa + 1, tempa + tempb + 1)) + list(range(1, tempa + 1))
                            + list(range(tempa + tempb + 1, tempa + tempb + tempx + 1))
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb + temgx, tuple(
                            list(range(temga + 1, temga + temgb + 1)) + list(range(1, temga + 1))
                            + list(range(temga + temgb + 1, temga + temgb + temgx + 1))
                        ))).evaluate()
                    if i.fields[1] < i.fields[2]:
                        tempx = i.fields[0].ind_s
                        tempa = i.fields[1].ind_s
                        tempb = i.fields[2].ind_s
                        temgx = i.fields[0].ind_g[0]
                        temga = i.fields[1].ind_g[0]
                        temgb = i.fields[2].ind_g[0]
                        i.fields = [i.fields[0], i.fields[2], i.fields[1]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + tempx, tuple(
                            list(range(1, tempx + 1))
                            + list(range(tempx + tempa + 1, tempx + tempa + tempb + 1))
                            + list(range(tempx + 1, tempx + tempa + 1))
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb + temgx, tuple(
                            list(range(1, temgx + 1))
                            + list(range(temgx + temga + 1, temgx + temga + temgb + 1))
                            + list(range(temgx + 1, temgx + temga + 1))
                        ))).evaluate()
                    if i.fields[0] < i.fields[2]:
                        tempa = i.fields[0].ind_s
                        tempx = i.fields[1].ind_s
                        tempb = i.fields[2].ind_s
                        temga = i.fields[0].ind_g[0]
                        temgx = i.fields[1].ind_g[0]
                        temgb = i.fields[2].ind_g[0]
                        i.fields = [i.fields[2], i.fields[0], i.fields[1]]
                        i.syms[0] = PermRule(i.syms[0], Permutation(tempa + tempb + tempx, tuple(
                            list(range(tempa + tempx + 1, tempa + tempx + tempb + 1))
                            + list(range(1, tempa + 1)) + list(range(tempa + 1, tempa + tempx + 1))
                        ))).evaluate()
                        i.syms[1] = PermRule(i.syms[1], Permutation(temga + temgb + temgx, tuple(
                            list(range(temga + temgx + 1, temga + temgx + temgb + 1))
                            + list(range(1, temga + 1)) + list(range(temga + 1, temga + temgx + 1))
                        ))).evaluate()
                    iter += 1
                    not_ordered = (i.fields[0] < i.fields[1] or i.fields[1] < i.fields[2]
                                   or i.fields[0] < i.fields[2])
                if iter >= 4:
                    raise Exception("IterationError")
                temp.append(i)
            else:
                raise Exception("ImplementationError")
        terms = temp

        """Symmetrize between commuting fields"""
        temp = []
        for i in terms:
            if len(i.fields) <= 1:
                temp.append(i)
            elif len(i.fields) == 2:
                if i.fields[0] == i.fields[1] and i.fields[1].kind == "c":
                    i.nmbr = i.nmbr.smul(Raw(1, 2, True))
                    t1 = copy.copy(i)
                    temp.append(t1)
                    t2 = copy.copy(i)
                    tempa = t2.fields[0].ind_s
                    tempb = t2.fields[1].ind_s
                    temga = t2.fields[0].ind_g[0]
                    temgb = t2.fields[1].ind_g[0]
                    t2.fields = [t2.fields[1], t2.fields[0]]
                    t2.syms[0] = PermRule(t2.syms[0], Permutation(tempa + tempb, tuple(
                        list(range(tempa + 1, tempa + tempb + 1)) + list(range(1, tempa + 1))
                    ))).evaluate()
                    t2.syms[1] = PermRule(t2.syms[1], Permutation(temga + temgb, tuple(
                        list(range(temga + 1, temga + temgb + 1)) + list(range(1, temga + 1))
                    ))).evaluate()
                    temp.append(t2)
                else:
                    temp.append(i)

            elif len(i.fields) == 3:
                if i.fields[0] == i.fields[1] == i.fields[2] and i.fields[1].kind == "c":
                    i.nmbr = i.nmbr.smul(Raw(1, 6, True))
                    t1 = copy.copy(i)
                    temp.append(t1)

                    t2 = copy.copy(i)
                    tempa = t2.fields[0].ind_s
                    tempb = t2.fields[1].ind_s
                    tempc = t2.fields[2].ind_s
                    temga = t2.fields[0].ind_g[0]
                    temgb = t2.fields[1].ind_g[0]
                    temgc = t2.fields[2].ind_g[0]
                    t2.fields = [t2.fields[1], t2.fields[2], t2.fields[0]]
                    t2.syms[0] = PermRule(t2.syms[0], Permutation(tempa + tempb + tempc, tuple([]
                        + list(range(tempa + 1, tempa + tempb + 1))
                        + list(range(tempa + tempb + 1, tempa + tempb + tempc + 1))
                        + list(range(1, tempa + 1))
                    ))).evaluate()
                    t2.syms[1] = PermRule(t2.syms[1], Permutation(temga + temgb + temgc, tuple([]
                        + list(range(temga + 1, temga + temgb + 1))
                        + list(range(temga + temgb + 1, temga + temgb + temgc + 1))
                        + list(range(1, temga + 1))
                    ))).evaluate()
                    temp.append(t2)

                    t3 = copy.copy(i)
                    tempa = t3.fields[0].ind_s
                    tempb = t3.fields[1].ind_s
                    tempc = t3.fields[2].ind_s
                    temga = t3.fields[0].ind_g[0]
                    temgb = t3.fields[1].ind_g[0]
                    temgc = t3.fields[2].ind_g[0]
                    t3.fields = [t3.fields[2], t3.fields[0], t3.fields[1]]
                    t3.syms[0] = PermRule(t3.syms[0], Permutation(tempa + tempb + tempc, tuple([]
                        + list(range(tempa + tempb + 1, tempa + tempb + tempc + 1))
                        + list(range(1, tempa + 1))
                        + list(range(tempa + 1, tempa + tempb + 1))
                    ))).evaluate()
                    t3.syms[1] = PermRule(t3.syms[1], Permutation(temga + temgb + temgc, tuple([]
                        + list(range(temga + temgb + 1, temga + temgb + temgc + 1))
                        + list(range(1, temga + 1))
                        + list(range(temga + 1, temga + temgb + 1))
                    ))).evaluate()
                    temp.append(t3)

                    t4 = copy.copy(i)
                    tempa = t4.fields[0].ind_s
                    tempb = t4.fields[1].ind_s
                    tempc = t4.fields[2].ind_s
                    temga = t4.fields[0].ind_g[0]
                    temgb = t4.fields[1].ind_g[0]
                    temgc = t4.fields[2].ind_g[0]
                    t4.fields = [t4.fields[0], t4.fields[2], t4.fields[1]]
                    t4.syms[0] = PermRule(t4.syms[0], Permutation(tempa + tempb + tempc, tuple([]
                        + list(range(1, tempa + 1))
                        + list(range(tempa + tempb + 1, tempa + tempb + tempc + 1))
                        + list(range(tempa + 1, tempa + tempb + 1))
                    ))).evaluate()
                    t4.syms[1] = PermRule(t4.syms[1], Permutation(temga + temgb + temgc, tuple([]
                        + list(range(1, temga + 1))
                        + list(range(temga + temgb + 1, temga + temgb + temgc + 1))
                        + list(range(temga + 1, temga + temgb + 1))
                    ))).evaluate()
                    temp.append(t4)

                    t5 = copy.copy(i)
                    tempa = t5.fields[0].ind_s
                    tempb = t5.fields[1].ind_s
                    tempc = t5.fields[2].ind_s
                    temga = t5.fields[0].ind_g[0]
                    temgb = t5.fields[1].ind_g[0]
                    temgc = t5.fields[2].ind_g[0]
                    t5.fields = [t5.fields[1], t5.fields[0], t5.fields[2]]
                    t5.syms[0] = PermRule(t5.syms[0], Permutation(tempa + tempb + tempc, tuple([]
                        + list(range(tempa + 1, tempa + tempb + 1))
                        + list(range(1, tempa + 1))
                        + list(range(tempa + tempb + 1, tempa + tempb + tempc + 1))
                    ))).evaluate()
                    t5.syms[1] = PermRule(t5.syms[1], Permutation(temga + temgb + temgc, tuple([]
                        + list(range(temga + 1, temga + temgb + 1))
                        + list(range(1, temga + 1))
                        + list(range(temga + temgb + 1, temga + temgb + temgc + 1))
                    ))).evaluate()
                    temp.append(t5)

                    t6 = copy.copy(i)
                    tempa = t6.fields[0].ind_s
                    tempb = t6.fields[1].ind_s
                    tempc = t6.fields[2].ind_s
                    temga = t6.fields[0].ind_g[0]
                    temgb = t6.fields[1].ind_g[0]
                    temgc = t6.fields[2].ind_g[0]
                    t6.fields = [t6.fields[2], t6.fields[1], t6.fields[0]]
                    t6.syms[0] = PermRule(t6.syms[0], Permutation(tempa + tempb + tempc, tuple([]
                        + list(range(tempa + tempb + 1, tempa + tempb + tempc + 1))
                        + list(range(tempa + 1, tempa + tempb + 1))
                        + list(range(1, tempa + 1))
                    ))).evaluate()
                    t6.syms[1] = PermRule(t6.syms[1], Permutation(temga + temgb + temgc, tuple([]
                        + list(range(temga + temgb + 1, temga + temgb + temgc + 1))
                        + list(range(temga + 1, temga + temgb + 1))
                        + list(range(1, temga + 1))
                    ))).evaluate()
                    temp.append(t6)

                elif i.fields[0] == i.fields[1] and i.fields[1].kind == "c":
                    i.nmbr = i.nmbr.smul(Raw(1, 2, True))
                    t1 = copy.copy(i)
                    temp.append(t1)

                    t5 = copy.copy(i)
                    tempa = t5.fields[0].ind_s
                    tempb = t5.fields[1].ind_s
                    tempc = t5.fields[2].ind_s
                    temga = t5.fields[0].ind_g[0]
                    temgb = t5.fields[1].ind_g[0]
                    temgc = t5.fields[2].ind_g[0]
                    t5.fields = [t5.fields[1], t5.fields[0], t5.fields[2]]
                    t5.syms[0] = PermRule(t5.syms[0], Permutation(tempa + tempb + tempc, tuple([]
                        + list(range(tempa + 1, tempa + tempb + 1))
                        + list(range(1, tempa + 1))
                        + list(range(tempa + tempb + 1, tempa + tempb + tempc + 1))
                    ))).evaluate()
                    t5.syms[1] = PermRule(t5.syms[1], Permutation(temga + temgb + temgc, tuple([]
                        + list(range(temga + 1, temga + temgb + 1))
                        + list(range(1, temga + 1))
                        + list(range(temga + temgb + 1, temga + temgb + temgc + 1))
                    ))).evaluate()
                    temp.append(t5)

                elif i.fields[1] == i.fields[2] and i.fields[1].kind == "c":
                    i.nmbr = i.nmbr.smul(Raw(1, 2, True))
                    t1 = copy.copy(i)
                    temp.append(t1)

                    t4 = copy.copy(i)
                    tempa = t4.fields[0].ind_s
                    tempb = t4.fields[1].ind_s
                    tempc = t4.fields[2].ind_s
                    temga = t4.fields[0].ind_g[0]
                    temgb = t4.fields[1].ind_g[0]
                    temgc = t4.fields[2].ind_g[0]
                    t4.fields = [t4.fields[0], t4.fields[2], t4.fields[1]]
                    t4.syms[0] = PermRule(t4.syms[0], Permutation(tempa + tempb + tempc, tuple([]
                        + list(range(1, tempa + 1))
                        + list(range(tempa + tempb + 1, tempa + tempb + tempc + 1))
                        + list(range(tempa + 1, tempa + tempb + 1))
                    ))).evaluate()
                    t4.syms[1] = PermRule(t4.syms[1], Permutation(temga + temgb + temgc, tuple([]
                        + list(range(1, temga + 1))
                        + list(range(temga + temgb + 1, temga + temgb + temgc + 1))
                        + list(range(temga + 1, temga + temgb + 1))
                    ))).evaluate()
                    temp.append(t4)
                else:
                    temp.append(i)
            else:
                raise Exception("UnknownError")
        terms = temp

        """Collect the coefficients to the number in front, again, just in case"""
        for i in terms:
            for j in i.consts:
                if not j.coef == Number(Raw(1, 1, True)):
                    i.nmbr = i.nmbr * j.coef
                    j.coef = Number(Raw(1, 1, True))
            for j in i.syms:
                if not j.coef == Number(Raw(1, 1, True)):
                    i.nmbr = i.nmbr * j.coef
                    j.coef = Number(Raw(1, 1, True))
            for j in i.fields:
                if not j.coef == Number(Raw(1, 1, True)):
                    i.nmbr = i.nmbr * j.coef
                    j.coef = Number(Raw(1, 1, True))

        """Bring the symbols to a linearly independent basis"""
        temp = []
        for i in terms:
            if i.syms[0] == catp:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(catn)
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(levi)
                temp2.nmbr = temp2.nmbr.smul(Raw(1,1,False))
                temp.append(temp2)
            elif i.syms[0] == catm:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(catn)
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(levi)
                temp2.nmbr = temp2.nmbr.smul(Raw(-1,1,False))
                temp.append(temp2)
            elif i.syms[0] == dogp:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(dogn)
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(levi)
                temp2.nmbr = temp2.nmbr.smul(Raw(1,1,False))
                temp.append(temp2)
            elif i.syms[0] == dogm:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(dogn)
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(levi)
                temp2.nmbr = temp2.nmbr.smul(Raw(-1,1,False))
                temp.append(temp2)
            elif i.syms[0] == foxp:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(foxn)
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(levi)
                temp2.nmbr = temp2.nmbr.smul(Raw(1,1,False))
                temp.append(temp2)
            elif i.syms[0] == foxm:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(foxn)
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(levi)
                temp2.nmbr = temp2.nmbr.smul(Raw(-1,1,False))
                temp.append(temp2)
            elif i.syms[0] == foxn:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(catn)
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(dogn)
                temp2.nmbr = temp2.nmbr.smul(Raw(-1,1,True))
                temp.append(temp2)
            elif i.syms[0] == owlp:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(eel)
                temp1.nmbr = temp1.nmbr.smul(Raw(-1,1,False))
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(elk)
                temp2.nmbr = temp2.nmbr.smul(Raw(1,1,True))
                temp.append(temp2)
            elif i.syms[0] == owlm:
                temp1 = copy.copy(i)
                temp1.syms[0] = copy.copy(eel)
                temp1.nmbr = temp1.nmbr.smul(Raw(-1,1,False))
                temp.append(temp1)
                temp2 = copy.copy(i)
                temp2.syms[0] = copy.copy(elk)
                temp2.nmbr = temp2.nmbr.smul(Raw(-1,1,True))
                temp.append(temp2)
            else:
                temp.append(i)
        terms = temp

        """Merge linearly dependent terms"""
        for ii,i in enumerate(terms):
            for jj,j in enumerate(terms):
                if jj <= ii:
                    pass
                else:
                    if i == j:
                        i.nmbr = i.nmbr + j.nmbr
                        j.nmbr = Number(Raw(0,1,True))
                    else:
                        pass
        temp = []
        for i in terms:
            if i.nmbr == Number(Raw(0,1,True)):
                pass
            else:
                temp.append(i)
        terms = temp

        return Expression(*terms, simp = simp0+1)

class Permutation:
    def __init__(self, num: int, elem: tuple):
        self.num = num
        self.elem = elem
        temp = list(elem[:])
        temp.sort()
        if not (len(elem) == num and temp == list(range(1, num+1))):
            raise Exception("PermutationError")

    def __str__(self):
        return str(self.elem)
    __repr__ = __str__

class PermRule:
    def __init__(self, symb: Symbol, perm: Permutation):
        self.symb = symb
        self.perm = perm

    def evaluate(self):
        if self.perm.elem == tuple(range(1, self.perm.num+1)):
            return copy.copy(self.symb)
        if self.symb.name in ["met", "delta"] and self.perm.elem == (2,1):
            return copy.copy(self.symb)

        if self.symb.name == "tiger":
            if self.perm.elem in [(1, 2, 3), (2, 3, 1), (3, 1, 2)]:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif self.perm.elem in [(1, 3, 2), (2, 1, 3), (3, 2, 1)]:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("RuleUnknown")

        if self.symb.name == "levi":
            if self.perm.elem in [(1, 2, 3, 4), (1, 3, 4, 2), (1, 4, 2, 3),
                                  (2, 1, 4, 3), (2, 3, 1, 4), (2, 4, 3, 1),
                                  (3, 1, 2, 4), (3, 2, 4, 1), (3, 4, 1, 2),
                                  (4, 1, 3, 2), (4, 2, 1, 3), (4, 3, 2, 1)]:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif self.perm.elem in [(1, 2, 4, 3), (1, 3, 2, 4), (1, 4, 3, 2),
                                    (2, 1, 3, 4), (2, 3, 4, 1), (2, 4, 1, 3),
                                    (3, 1, 4, 2), (3, 2, 1, 4), (3, 4, 2, 1),
                                    (4, 1, 2, 3), (4, 2, 3, 1), (4, 3, 1, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("RuleUnknown")

        if self.symb.name == "catn":
            if self.perm.elem in [(1, 2, 3, 4), (2, 1, 4, 3), (3, 4, 1, 2), (4, 3, 2, 1)]:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif self.perm.elem in [(1, 2, 4, 3), (2, 1, 3, 4), (3, 4, 2, 1), (4, 3, 1, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif self.perm.elem in [(1, 4, 3, 2), (2, 3, 4, 1), (3, 2, 1, 4), (4, 1, 2, 3)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogn"
                return tcopy
            elif self.perm.elem in [(1, 4, 2, 3), (2, 3, 1, 4), (3, 2, 4, 1), (4, 1, 3, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogn"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif self.perm.elem in [(1, 3, 2, 4), (2, 4, 1, 3), (3, 1, 4, 2), (4, 2, 3, 1)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxn"
                return tcopy
            elif self.perm.elem in [(1, 3, 4, 2), (2, 4, 3, 1), (3, 1, 2, 4), (4, 2, 1, 3)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxn"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("RuleUnknown")

        if self.symb.name == "dogn":
            if self.perm.elem in [(1, 2, 3, 4), (2, 1, 4, 3), (3, 4, 1, 2), (4, 3, 2, 1)]:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif self.perm.elem in [(1, 3, 2, 4), (2, 4, 1, 3), (3, 1, 4, 2), (4, 2, 3, 1)]:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif self.perm.elem in [(1, 4, 3, 2), (2, 3, 4, 1), (3, 2, 1, 4), (4, 1, 2, 3)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catn"
                return tcopy
            elif self.perm.elem in [(1, 3, 4, 2), (2, 4, 3, 1), (3, 1, 2, 4), (4, 2, 1, 3)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catn"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif self.perm.elem in [(1, 4, 2, 3), (2, 3, 1, 4), (3, 2, 4, 1), (4, 1, 3, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxn"
                return tcopy
            elif self.perm.elem in [(1, 2, 4, 3), (2, 1, 3, 4), (3, 4, 2, 1), (4, 3, 1, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxn"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("RuleUnknown")

        if self.symb.name == "foxn":
            if self.perm.elem in [(1, 2, 3, 4), (2, 1, 4, 3), (3, 4, 1, 2), (4, 3, 2, 1)]:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif self.perm.elem in [(1, 4, 3, 2), (2, 3, 4, 1), (3, 2, 1, 4), (4, 1, 2, 3)]:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif self.perm.elem in [(1, 3, 2, 4), (2, 4, 1, 3), (3, 1, 4, 2), (4, 2, 3, 1)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catn"
                return tcopy
            elif self.perm.elem in [(1, 4, 2, 3), (2, 3, 1, 4), (3, 2, 4, 1), (4, 1, 3, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catn"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif self.perm.elem in [(1, 3, 4, 2), (2, 4, 3, 1), (3, 1, 2, 4), (4, 2, 1, 3)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogn"
                return tcopy
            elif self.perm.elem in [(1, 2, 4, 3), (2, 1, 3, 4), (3, 4, 2, 1), (4, 3, 1, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogn"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("RuleUnknown")

        if self.symb.name == "catp":
            temp1 = PermRule(catn, self.perm).evaluate().name
            temp2 = PermRule(catn, self.perm).evaluate().coef == Number(Raw(1,1,True))
            temp3 = PermRule(levi, self.perm).evaluate().coef == Number(Raw(1,1,True))
            if temp1 == "catn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif temp1 == "dogn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                return tcopy
            elif temp1 == "foxn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                return tcopy
            elif temp1 == "catn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                return tcopy
            elif temp1 == "dogn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                return tcopy
            elif temp1 == "foxn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                return tcopy
            elif temp1 == "catn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("UnknownError")

        if self.symb.name == "catm":
            temp1 = PermRule(catn, self.perm).evaluate().name
            temp2 = PermRule(catn, self.perm).evaluate().coef == Number(Raw(1,1,True))
            temp3 = PermRule(levi, self.perm).evaluate().coef == Number(Raw(1,1,True))
            if temp1 == "catn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif temp1 == "dogn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                return tcopy
            elif temp1 == "foxn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                return tcopy
            elif temp1 == "catn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                return tcopy
            elif temp1 == "dogn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                return tcopy
            elif temp1 == "foxn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                return tcopy
            elif temp1 == "catn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("UnknownError")

        if self.symb.name == "dogp":
            temp1 = PermRule(dogn, self.perm).evaluate().name
            temp2 = PermRule(dogn, self.perm).evaluate().coef == Number(Raw(1,1,True))
            temp3 = PermRule(levi, self.perm).evaluate().coef == Number(Raw(1,1,True))
            if temp1 == "dogn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif temp1 == "catn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                return tcopy
            elif temp1 == "foxn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                return tcopy
            elif temp1 == "dogn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                return tcopy
            elif temp1 == "catn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                return tcopy
            elif temp1 == "foxn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                return tcopy
            elif temp1 == "dogn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("UnknownError")

        if self.symb.name == "dogm":
            temp1 = PermRule(dogn, self.perm).evaluate().name
            temp2 = PermRule(dogn, self.perm).evaluate().coef == Number(Raw(1,1,True))
            temp3 = PermRule(levi, self.perm).evaluate().coef == Number(Raw(1,1,True))
            if temp1 == "dogn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif temp1 == "catn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                return tcopy
            elif temp1 == "foxn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                return tcopy
            elif temp1 == "dogn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                return tcopy
            elif temp1 == "catn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                return tcopy
            elif temp1 == "foxn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                return tcopy
            elif temp1 == "dogn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("UnknownError")

        if self.symb.name == "foxp":
            temp1 = PermRule(foxn, self.perm).evaluate().name
            temp2 = PermRule(foxn, self.perm).evaluate().coef == Number(Raw(1,1,True))
            temp3 = PermRule(levi, self.perm).evaluate().coef == Number(Raw(1,1,True))
            if temp1 == "foxn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif temp1 == "catn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                return tcopy
            elif temp1 == "dogn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                return tcopy
            elif temp1 == "foxn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                return tcopy
            elif temp1 == "catn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                return tcopy
            elif temp1 == "dogn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                return tcopy
            elif temp1 == "foxn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("UnknownError")

        if self.symb.name == "foxm":
            temp1 = PermRule(foxn, self.perm).evaluate().name
            temp2 = PermRule(foxn, self.perm).evaluate().coef == Number(Raw(1,1,True))
            temp3 = PermRule(levi, self.perm).evaluate().coef == Number(Raw(1,1,True))
            if temp1 == "foxn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif temp1 == "catn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                return tcopy
            elif temp1 == "dogn" and temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                return tcopy
            elif temp1 == "foxn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                return tcopy
            elif temp1 == "catn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                return tcopy
            elif temp1 == "dogn" and temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                return tcopy
            elif temp1 == "foxn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "foxp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogp"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "foxn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "catn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "catm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif temp1 == "dogn" and not temp2 and not temp3:
                tcopy = copy.copy(self.symb)
                tcopy.name = "dogm"
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("UnknownError")

        if self.symb.name in ["eel", "elk", "owlp", "owlm"]:
            if self.perm.elem in [(1, 2, 3, 4, 5, 6), (3, 4, 5, 6, 1, 2), (5, 6, 1, 2, 3, 4)]:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif self.perm.elem in [(1, 2, 5, 6, 3, 4), (3, 4, 1, 2, 5, 6), (5, 6, 3, 4, 1, 2)]:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            elif self.perm.elem in [(1, 2, 3, 4, 5, 6), (1, 2, 4, 3, 6, 5),
                                    (2, 1, 3, 4, 6, 5), (2, 1, 4, 3, 5, 6)]:
                tcopy = copy.copy(self.symb)
                return tcopy
            elif self.perm.elem in [(2, 1, 3, 4, 5, 6), (1, 2, 4, 3, 5, 6),
                                    (1, 2, 3, 4, 6, 5), (2, 1, 4, 3, 6, 5)]:
                tcopy = copy.copy(self.symb)
                tcopy.coef = tcopy.coef.smul(-1)
                return tcopy
            else:
                raise Exception("RuleUnknown")

        raise Exception("RuleUnknown")

dermat = Field("P", kind="d", ind_s=1, ind_g=(0,))

class Delphi:

    def __init__(self):
        self.rec = Listener({})
        Speaker().intro()
        self.main_menu()

    def main_menu(self):
        Speaker().main_menu()
        entry = Speaker().inputter(lambda s: s in ["a", "b", "c"])
        if entry == "a":
            self.guide()
        elif entry == "b":
            self.about()
        else:
            pass

    def about(self):
        Speaker().about()
        Speaker().back_to_main()
        self.main_menu()

    def guide(self):
        self.fragments()
        self.background()
        self.gauge()
        self.freezing()
        self.affine()
        Speaker().before_action()
        self.organizer()

    def fragments(self):
        Speaker().fragments()
        entry = Speaker().inputter(lambda s: s in ["1", "2", "3", "4", "5", "6", "7", "8", "9"])
        entry = int(entry)
        self.rec.dct["frag_num"] = entry

    def background(self):
        Speaker().background(self.rec.dct["frag_num"])
        def temp_func(s):
            if s == "0":
                return True
            temp = s.split(",")
            try:
                temp = [int(i) for i in temp]
                temp = all([i >= 1 and i <= self.rec.dct["frag_num"] for i in temp])
                return temp
            except ValueError:
                return False
        entry = Speaker().inputter(temp_func)
        if entry == "0":
            self.rec.dct["frag_der"] = ()
        else:
            entry = entry.split(",")
            entry = list(set([int(i) for i in entry]))
            entry.sort()
            entry = tuple(entry)
            self.rec.dct["frag_der"] = entry

    def gauge(self):
        Speaker().gauge()

    def freezing(self):
        Speaker().freezing()
        newf = Speaker().inputter(lambda s: s in ["y", "n"])
        self.rec.dct["fields"] = []
        while newf == "y":
            self.field_designer()
            Speaker().field_new()
            newf = Speaker().inputter(lambda s: s in ["y", "n"])

    def field_designer(self):
        Speaker().field_design_1()
        name = Speaker().namer()
        Speaker().field_design_2()
        lor = Speaker().inputter(lambda s: s in ["s", "v", "t"])
        if lor == "s":
            Speaker().field_design_2_s()
            chir = Speaker().inputter(lambda s: s in ["0", "p", "l", "r"])
        elif lor == "v":
            Speaker().field_design_2_v()
            chir = Speaker().inputter(lambda s: s in ["0", "p", "l", "r"])
        else:
            Speaker().field_design_2_t()
            chir = Speaker().inputter(lambda s: s in ["0", "p", "l", "r"])
        Speaker().field_design_3()
        rep = Speaker().inputter(lambda s: s in ["t", "a"])
        Speaker().field_design_4(self.rec.dct["frag_num"])
        degen = Speaker().num_lister(self.rec.dct["frag_num"])
        Speaker().field_design_5a()
        Displayer().field_displayer(name, lor, chir, rep, degen)
        Speaker().field_design_5b()
        agree = Speaker().inputter(lambda s: s in ["y", "n"])
        if agree == "y":
            self.rec.dct["fields"] = self.rec.dct["fields"] + [{
                "name": name, "lor": lor, "chir": chir, "rep": rep, "degen": degen
            }]
        else:
            return None

    def affine(self):
        Speaker().affine()

    def organizer(self):
        """self.rec.dct["frag_num"] = 1
        self.rec.dct["frag_der"] = (1,)
        self.rec.dct["fields"] = [
            {"name": "A", "lor": "v", "chir": "0", "rep": "a", "degen": [Number(Raw(1,1,True))]},
            {"name": "B", "lor": "t", "chir": "p", "rep": "a", "degen": [Number(Raw(1,1,True))]}
        ]"""

        frag_num = self.rec.dct["frag_num"]
        frag_der = self.rec.dct["frag_der"]
        fields = self.rec.dct["fields"]
        dct_ind_s = {"s": 0, "v": 1, "t": 2}
        dct_ind_g = {"t": 0, "a": 1}
        dct_Gamma = {"sl": Gamma(c1=Number(Raw(1,1,True))),
                     "sr": Gamma(c2=Number(Raw(1,1,True))),
                     "vl": Gamma(c3=Number(Raw(1,1,True))),
                     "vr": Gamma(c4=Number(Raw(1,1,True))),
                     "tl": Gamma(c5=Number(Raw(1,1,True))),
                     "tr": Gamma(c6=Number(Raw(1,1,True)))}
        dct_Gener = {"t": Gener(c1=Number(Raw(1,1,True))),
                     "a": Gener(c2=Number(Raw(1,1,True)))}
        terms = []
        for i in range(frag_num):
            list_at_frag = []
            if i+1 in frag_der:
                list_at_frag.append([["vr", "t"], dermat])
                list_at_frag.append([["vl", "t"], dermat])
            for j in fields:
                if not (j["degen"][i] == Number(Raw(0, 1, True))):
                    if j["chir"] in ["l", "r"]:
                        list_at_frag.append(
                            [[j["lor"] + j["chir"], j["rep"]],
                             Field(j["name"], kind="c", der=0, coef = j["degen"][i],
                                   ind_s=dct_ind_s[j["lor"]], ind_sp=0,
                                   ind_g=(dct_ind_g[j["rep"]],), ind_gp=0)]
                        )
                    elif j["chir"] == "0":
                        list_at_frag.append(
                            [[j["lor"] + "r", j["rep"]],
                             Field(j["name"], kind="c", der=0, coef = j["degen"][i],
                                   ind_s=dct_ind_s[j["lor"]], ind_sp=0,
                                   ind_g=(dct_ind_g[j["rep"]],), ind_gp=0)]
                        )
                        list_at_frag.append(
                            [[j["lor"] + "l", j["rep"]],
                             Field(j["name"], kind="c", der=0, coef = j["degen"][i],
                                   ind_s=dct_ind_s[j["lor"]], ind_sp=0,
                                   ind_g=(dct_ind_g[j["rep"]],), ind_gp=0)]
                        )
                    elif j["chir"] == "p":
                        list_at_frag.append(
                            [[j["lor"] + "r", j["rep"]],
                             Field(j["name"], kind="c", der=0, coef = j["degen"][i],
                                   ind_s=dct_ind_s[j["lor"]], ind_sp=0,
                                   ind_g=(dct_ind_g[j["rep"]],), ind_gp=0)]
                        )
                        list_at_frag.append(
                            [[j["lor"] + "l", j["rep"]],
                             Field(j["name"], kind="c", der=0, coef = j["degen"][i].smul(-1),
                                   ind_s=dct_ind_s[j["lor"]], ind_sp=0,
                                   ind_g=(dct_ind_g[j["rep"]],), ind_gp=0)]
                        )
                    else:
                        raise Exception("LogicalError")
            triplets_at_frag = [[[dct_Gamma[i[0][0]], dct_Gamma[j[0][0]], dct_Gamma[k[0][0]]],
                                 [dct_Gener[i[0][1]], dct_Gener[j[0][1]], dct_Gener[k[0][1]]],
                                 [i[1], j[1], k[1]]]
                                for i in list_at_frag for j in list_at_frag for k in list_at_frag]
            for i in triplets_at_frag:
                i[0] = i[0][0].trimul(i[0][1], i[0][2])
                i[1] = i[1][0].trimul(i[1][1], i[1][2])
            temp = []
            for i in triplets_at_frag:
                for j in i[0]:
                    for k in i[1]:
                        temp.append([j, k, i[2]])
            triplets_at_frag = temp
            terms_at_frag = [Term(i[0][0] * i[1][0], [], [i[0][1], i[1][1]], i[2])
                             for i in triplets_at_frag]
            terms += terms_at_frag
        answer = Expression(*terms)
        print(answer)
        return answer

class Environment:
    def __init__(self, dim=4):
        self.dimension = dim
        if dim == 4:
            self.levi = "bilmemne"
        pass

class Listener:
    def __init__(self, dct: dict):
        self.dct = dct

class Speaker:
    def __init__(self):
        pass

    @staticmethod
    def inputter(cond):
        temp = "unregistered"
        while not cond(temp):
            temp = input()
        return temp

    @staticmethod
    def namer():
        name = ""
        ticket = False
        while ticket == False:
            name = input()
            forbidden = ["", "P", "i", "f", "h", "\u03B3", "\u03B4", "\u03B5", "\u03B7",
                         "one", "met", "lev", "levi", "cat", "catn", "catp", "catm",
                         "dog", "dogn", "dogp", "dogm", "fox", "foxn", "foxp", "foxm",
                         "eel", "elk", "owl", "owlp", "owlm",
                         "delta", "tiger", "hyena"]
            forbidden_chars = [" ", "_", ",", "(", ")", "/", "\u2202"]
            if name == "":
                print("The name cannot be an empty string.")
            elif len(name) > 9:
                print("This name is too long. Please try again.")
            elif not name.isalpha():
                print("The name has to be alphabetical.")
            elif not all([(i not in name) for i in forbidden_chars]):
                print("This name contains a reserved character. Please try again.")
            elif name in forbidden:
                print("This name is reserved in our code for another purpose. Please try again.")
            else:
                ticket = True
        return name

    @staticmethod
    def num_lister(n: int):
        answer = []
        ticket = False
        while ticket == False:
            entry = input()
            if entry == "":
                continue
            entry = entry.split(",")
            if len(entry) != n:
                continue
            entry = [i.split("/") for i in entry]
            if not all([len(i) <= 2 for i in entry]):
                continue
            if not all([all([len(j) > 0 for j in i]) for i in entry]):
                continue
            entry = [i + ["1"] if len(i) == 1 else i for i in entry]
            entry = [i + [True] for i in entry]
            entry = [[i[0][:-1], i[1], False] if i[0][-1] == "i" else i for i in entry]
            if not all([len(i[0]) > 0 for i in entry]):
                continue
            try:
                entry = [[int(i[0]), int(i[1]), i[2]] for i in entry]
            except ValueError:
                continue
            if not all([i[1] != 0 for i in entry]):
                print("Zero division error. Please try again.")
                continue
            answer = [Number(Raw(*i)) for i in entry]
            ticket = True
        return answer

    @staticmethod
    def intro():
        print(" ")
        print("Welcome to Delphi, v0.1 @ Microsoft")
        print(" ")
        print("Press Enter to continue...")
        input()
        print(" ")
        print("This program creates theories for fundamental physics " +
              "from freezings of the cubic matrix action S = \u2159 Tr(M\u00B3).")
        print(" ")
        print("While the overall shape of the action is fixed, "
              + "the set of degrees of freedom inside it is not.")
        print("Various parts of the matrix M can be chosen to be variable or fixed, "
              + "which we call a 'freezing'.")
        print("Different freezing selections gives rise to different theories.")
        print(" ")
        print("You can create alternative universes here, "
              + "or search for an Ultimate Theory of Everything that describes yours.")
        print(" ")
        input()

    @staticmethod
    def main_menu():
        print("Main Menu")
        print("[a] Start")
        print("[b] About")
        print("[c] Quit")
        print(" ")

    @staticmethod
    def back_to_main():
        print("Press Enter to go back to the main menu.")
        input()

    @staticmethod
    def about():
        print("To be implemented...")
        print(" ")

    @staticmethod
    def fragments():
        print(" ")
        print("Part 1 of 3: The overall structure of M.")
        print(" ")
        input()
        print("In most cases, the matrix M is block-diagonal. We call each block a 'fragment'.")
        print(" ")
        print("The fragments decouple in the action, but a single degree of freedom can stretch")
        print("across multiple fragments, which adds a 'degeneracy' to the freezing options.")
        print(" ")
        print("How many fragments should M consist of? [1..9]")

    @staticmethod
    def background(fr_n):
        print(" ")
        print("The matrix M is required to be invariant under the Lorentz group action of SO(3,1),")
        print("which enables a 3+1-dimensional spacetime interpretation for the resulting actions.")
        print("We prepared this background structure for your theory.")
        input()
        print("Since the Lorentz group is reducible, most invariant orbits have a further 'degeneracy'.")
        print("The 6-dim Lorentz eigenspace can be considered in three parts:")
        print("2-dim space of scalars [s0], pseudo-scalars [sp], "
              + "left-handed scalars [sl], right-handed scalars [sr].")
        print("2-dim space of vectors [v0], pseudo-vectors [vp], "
              + "left-handed vectors [vl], right-handed vectors [vr].")
        print("2-dim space of skew-symmetric 2-tensors, which can be non-chiral [t0], inverted [tp], "
              + "self-dual [tl], or anti-self-dual [tr].")
        input()
        print("The derivatives in the action originate from certain blocks in the matrix M, "
              + "whose values are fixed to the")
        print("so-called 'derivative matrices'. The derivative matrices are diagonal and "
              + "their spectrum is the set of integers.")
        input()
        print(f"Which of the {fr_n} fragments should contain these derivative matrices?")
        print(f"(Enter digits from 1 to {fr_n} separated by commas. Enter 0 for no derivatives.)")

    @staticmethod
    def gauge():
        print(" ")
        print("The matrix M is also required to be invariant under the group action of a set of gauge groups.")
        input()
        print("In the current version of Delphi, the options are limited to a single SU(n) group with n>1,")
        print("and the trivial [t] and adjoint [a] representations.")
        input()

    @staticmethod
    def freezing():
        print(" ")
        print("Part 2 of 3: Freezing.")
        print(" ")
        input()
        print("The matrix is built, now you can choose what parts of it should be variable in the action.")
        input()
        print("The selected parts must be invariant under the symmetry group actions on M. Therefore the options")
        print("are labeled by their representations. Each such choice will become a distinct 'field' in the theory.")
        input()
        print("Would you like to add a new field? [y/n]")

    @staticmethod
    def field_new():
        print(" ")
        print("Would you like to add a new field? [y/n]")

    @staticmethod
    def field_design_1():
        print("How should this field be named? This determines how it will be displayed in the final action.")

    @staticmethod
    def field_design_2():
        print("Is it a scalar, a vector, or a 2-form? [s/v/t]")

    @staticmethod
    def field_design_2_s():
        print("Is it a true scalar [0], a pseudo-scalar [p], a left-handed [l], or a right-handed combination [r]?")

    @staticmethod
    def field_design_2_v():
        print("Is it a true vector [0], a pseudo-vector [p], a left-handed vector [l], or a right-handed vector [r]?")

    @staticmethod
    def field_design_2_t():
        print("Is it non-chiral [0], Hodge-inverted non-chiral [p], self-dual [l], or anti-self-dual [r]?")

    @staticmethod
    def field_design_3():
        print("Is it in the trivial or adjoint representation of SU(n)? [t/a]")

    @staticmethod
    def field_design_4(fr_n):
        print("Please specify the degeneracy of this field across the fragments.")
        print(f"The entry needs to be {fr_n} integers or integer fractions ('m/n') separated by commas.")

    @staticmethod
    def field_design_5a():
        print("You selected to unfreeze:")

    @staticmethod
    def field_design_5b():
        print("Do you agree? [y/n]")

    @staticmethod
    def affine():
        print(" ")
        print("Part 3 of 3: Affine parameters.")
        print(" ")
        print("This is yet to be implemented.")
        input()

    @staticmethod
    def before_action():
        print(" ")
        print("The selection is complete. Press Enter to see the action.")
        print(" ")
        input()

class Displayer:
    def __init__(self):
        pass

    unicode_subscript_digit = {
        0: "\u2080", 1: "\u2081", 2: "\u2082", 3: "\u2083", 4: "\u2084",
        5: "\u2085", 6: "\u2086", 7: "\u2087", 8: "\u2088", 9: "\u2089"
    }

    @staticmethod
    def field_displayer(name, lor, chir, rep, degen):
        answer = ""
        d = [(i, j) for i, j in enumerate(degen) if not (j == Number(Raw(0, 1, True)))]
        if len(d) == 0:
            print(0)
            return None
        elif len(d) == 1:
            dd = d[0]
            temp = str(dd[1]) + " " + "\u018E" + Displayer().unicode_subscript_digit[dd[0]+1]
            answer += temp
        else:
            temp = ""
            for dd in d:
                temp += " + " + str(dd[1]) + " \u018E" + Displayer().unicode_subscript_digit[dd[0]+1]
            temp = "(" + temp[3:] + ")"
            answer += temp
        answer += " \u2297 "
        if lor == "s" and chir == "0":
            answer += "1" + Displayer().unicode_subscript_digit[4]
        elif lor == "s" and chir == "p":
            answer += "\u03B3" + Displayer().unicode_subscript_digit[5]
        elif lor == "s" and chir == "l":
            answer += "\u03B3\u2097"
        elif lor == "s" and chir == "r":
            answer += "\u03B3\u1d63"
        elif lor == "v" and chir == "0":
            answer += "\u03B3\u1d45"
        elif lor == "v" and chir == "p":
            answer += "\u03B3\u1d45" + "\u03B3" + Displayer().unicode_subscript_digit[5]
        elif lor == "v" and chir == "l":
            answer += "\u03B3\u1d45" + "\u03B3\u2097"
        elif lor == "v" and chir == "r":
            answer += "\u03B3\u1d45" + "\u03B3\u1d63"
        elif lor == "t" and chir == "0":
            answer += "\u0393\u1d45\u1d5d"
        elif lor == "t" and chir == "p":
            answer += "\u0393\u1d45\u1d5d" + "\u03B3" + Displayer().unicode_subscript_digit[5]
        elif lor == "t" and chir == "l":
            answer += "\u0393\u1d45\u1d5d" + "\u03B3\u2097"
        elif lor == "t" and chir == "r":
            answer += "\u0393\u1d45\u1d5d" + "\u03B3\u1d63"
        else:
            raise Exception("UnknownError")
        answer += " \u2297 "
        if rep == "t":
            answer += "1"
        elif rep == "a":
            answer += "\u03C4\u02b2"
        else:
            raise Exception("UnknownError")
        answer += " \u2297 "
        answer += name
        if lor == "v":
            answer += "\u1d45"
        elif lor == "t":
            answer += "\u1d45\u1d5d"
        if rep == "a":
            answer += "\u02b2"
        answer += " \u2208 M"
        print(answer)
        return None


Oracle = Delphi()
#Oracle.organizer()
