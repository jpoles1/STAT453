#Contingency Table Module - Provides an object representing the contigency table and many related analytic functions
#Written by Jordan Poles
#Provided under the MIT License (c) 2017 Jordan Poles 
from IPython.display import display
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from re import search
from math import exp, sqrt, floor, log

class ContTable:
    def __init__(self, de, nde, dne, ndne, verbose=1, xlabs=["D", "~D"], ylabs=["E", "~E"]):
        assert de >= 0 and nde >= 0 and dne >= 0 and ndne >= 0, "All values in the contingency table must be positive."
        self.desc = [
            "Diseased and exposed", 
            "Not disease but exposed", 
            "Diseased but not exposed", 
            "Not diseased and not exposed"
        ]
        self.xlabels = xlabs
        self.ylabels = ylabs
        self.de = de; #Diseased and exposed = a
        self.a = de;
        self.nde = nde; #Not disease but exposed = b
        self.b = nde;
        self.dne = dne; #Diseased but not exposed = c
        self.c = dne;
        self.ndne = ndne; #Not diseased and not exposed = d
        self.d = ndne;
        self.num_d = de + dne;
        self.num_nd = nde + ndne;
        self.num_e = de + nde;
        self.num_ne = dne + ndne;
        self.num_tot = de+nde+dne+ndne
        self.cont_matrix = np.matrix([[self.de,self.nde],[self.dne,self.ndne]])
        self.cont_table = pd.DataFrame(self.cont_matrix)
        self.cont_table.columns = self.ylabels
        self.cont_table.index = self.xlabels
        if verbose:
            self.display()
            self.summarize(verbose=verbose)
    def relative_risk(self, verbose=1):
        #Basis for multiplicative risk
        #RR = P(D|E)/P(D|~E)
        num = self.de/self.num_e
        denom = self.dne/self.num_ne
        self.RR = num/denom
        if verbose:
            print("Relative Risk (RR):", self.RR)
        if verbose > 1:
            if self.RR == 1:
                print("Relative Risk (RR) indicates independence between D & E.")
            elif self.RR > 1:
                print("=> Relative Risk (RR) indicates a positive relationship between D & E.")
            else:
                print("=> Relative Risk (RR) indicates a negative relationship between D & E.")
            print("")
        assert self.RR >= 0, "Relative risk must be non-negative!"
        #Restriction on range (below) becomes problematic with common disease outcomes
        assert self.RR <= 1/denom, "Relative Risk is always <= 1/P(D|~E)"
        return self.RR
    def odds_ratio(self, CI=.95, verbose=1):
        #OR = [P(D|E)/P(~D|E)]/[P(D|~E)/P(~D|~E)] = (ad)/(bc)
        self.OR = (self.de*self.ndne)/(self.nde*self.dne);
        self.logOR = log(self.OR)
        z = -norm.interval(CI)[0]
        logORVar = (1/self.a) + (1/self.b) + (1/self.c) + (1/self.d)
        #Upper bound
        upper = exp(self.logOR + z * sqrt(logORVar))
        #Lower bound
        lower = exp(self.logOR - z * sqrt(logORVar))
        if verbose:
            print("Odds Ratio (OR):", self.OR)
            if CI:
                print("=> {:d}% CI using logOR: {:f}".format(floor(CI*100), self.logOR)) 
                print("=> => logOR Var:", logORVar)
                print("=> => Upper Bound:", upper)
                print("=> => Lower Bound:", lower)
        if verbose > 1:
            if self.OR == 1:
                print("=> Odds Ratio (OR) indicates independence between D & E.")
            elif self.OR > 1:
                print("=> Odds Ratio (OR) indicates a positive relationship between D & E.")
            else:
                print("=> Odds Ratio (OR) indicates a negative relationship between D & E.")
            print("")
        assert self.OR >= 0, "Odds Ratio must be non-negative!"
        #Odds Ratio has no upper bound, unlike Relative Risk
        return dict(OR=self.OR, logOR=self.logOR, bounds=(lower, upper))
    def ss_odds_ratio(self, CI=.95, verbose=1):
        self.ss_OR = (self.a*self.d)/((self.b + 1)*(self.c + 1))
        a = self.a + .5
        b = self.b + .5
        c = self.c + .5
        d = self.d + .5
        self.ss_logOR = log((a*d)/(b*c))
        z = -norm.interval(CI)[0]
        ss_logORVar = (1/a) + (1/b) + (1/c) + (1/d)
        #Upper bound
        upper = exp(self.ss_logOR + z * sqrt(ss_logORVar))
        #Lower bound
        lower = exp(self.ss_logOR - z * sqrt(ss_logORVar))
        if verbose:
            print("Small Sample Odds Ratio (ss_OR):", self.ss_OR)
            if CI:
                print("=> {:d}% CI using ss_logOR: {:f}".format(floor(CI*100), self.ss_logOR))
                print("=> => ss_logOR Var:", ss_logORVar)
                print("=> => Upper Bound:", upper)
                print("=> => Lower Bound:", lower)
        if verbose > 1:
            if self.OR == 1:
                print("=> Small Sample Odds Ratio (OR) indicates independence between D & E.")
            elif self.OR > 1:
                print("=> Small Sample Odds Ratio (OR) indicates a positive relationship between D & E.")
            else:
                print("=> Small Sample Odds Ratio (OR) indicates a negative relationship between D & E.")
            print("")
        return dict(OR=self.OR, logOR=self.logOR, bounds=(lower, upper))
    def excess_risk(self, verbose=1):
        #Basis for additive risk
        #ER = P(D|E) - P(D|~E) 
        first = self.de/self.num_e
        second = self.dne/self.num_ne
        self.ER = first-second
        if verbose:
            print("Excess Risk (ER):", self.ER)
        if verbose > 1:
            if self.ER == 0:
                print("=> Excess Risk (ER) indicates independence between D & E.")
            elif self.ER > 0:
                print("=> Excess Risk (ER) indicates a positive relationship between D & E.")
            else:
                print("=> Excess Risk (ER) indicates a negative relationship between D & E.")
            print("")
        assert self.ER >= -1 and self.ER <= 1, "Excess Risk (ER) is bounded by [-1, 1]." 
        return self.ER
    def attributable_risk(self, verbose=1):
        #AR = [P(D) - P(D|~E)]/P(D)
        pd = (self.de + self.dne)/self.num_tot
        pdne = self.dne/self.num_ne
        self.AR = (pd - pdne)/pd
        if verbose:
            print("Attributable Risk (AR):", self.AR)
        if verbose > 1:
            if self.AR == 0:
                print("=> Attributable Risk (AR) indicates independence between D & E.")
            elif self.AR > 0:
                print("=> Attributable Risk (AR) indicates a positive relationship between D & E.")
            else:
                print("=> Attributable Risk (AR) indicates a negative relationship between D & E.")
            print("")
        return self.AR
    def chisq_indep(self, alpha=.05, verbose=1):
        #Statistic = n(ad-bc)^2/(a+b)(a+c)(b+d)(c+d)
        #a=de; b=nde; c=dne; d=ndne
        #Large values suggest that the factors are not indep.
        num = (self.de*self.ndne)-(self.nde*self.dne)
        var_hat = (self.de+self.nde)*(self.de+self.dne)*(self.nde+self.ndne)*(self.dne+self.ndne)/self.num_tot
        test_statistic = np.square(num)/var_hat
        pval = 1-chi2.cdf(test_statistic, 1)
        if verbose:
            print("ChiSquared Independence Test (n=%i)"%self.num_tot)
            print("=> test statistic =", test_statistic)
            print("=> p-value ≈", pval)
        if verbose > 1:
            if pval < alpha:
                print("=> ChiSquared Test suggests an association between D & E.")
            else:
                print("=> ChiSquared Test indicates independence between D & E.")
            print("")
        return pval
    def display(self):
        print("2x2 Contingency Table:")
        display(self.cont_table)
    def summarize(self, verbose=1):
        self.relative_risk(verbose=verbose)
        self.odds_ratio(verbose=verbose)
        self.excess_risk(verbose=verbose)
        self.attributable_risk(verbose=verbose)
        self.chisq_indep(verbose=verbose)
    #Utility Function
    #Takes a string of form: "P(__)" where the blank may be filled with D, D|E, E|D, ~E, etc.
    #Eg: self.get_prob("P(~E|D)")
    def get_prob(self, prob_string):
        valid = search("P\((.*)\)", prob_string)
        if valid:
            cond_prob = search("([\~a-zA-Z]*)\|([\~a-zA-Z]*)", valid.group(1))
            if cond_prob:
                a = cond_prob.group(1)
                b = cond_prob.group(2)
                if a in self.xlabels:
                    return self.cont_table.loc[a, b]/self.cont_table.loc[:, b].sum()
                else:
                    return self.cont_table.loc[b, a]/self.cont_table.loc[:, a].sum()
            else:
                prob_chars = valid.group(1)
                if prob_chars in self.xlabels:
                    return self.cont_table.loc[prob_chars, :].sum()/self.num_tot
                else:
                    return self.cont_table.loc[:, prob_chars].sum()/self.num_tot