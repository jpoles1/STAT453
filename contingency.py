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
        self.relative_risk(verbose=verbose)
        self.odds_ratio(verbose=verbose)
        self.excess_risk(verbose=verbose)
        self.attributable_risk(verbose=verbose)
        self.chisq_indep(verbose=verbose)
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
            print("=> p-value =", pval)
        if verbose > 1:
            if pval < .05:
                print("=> ChiSquared Test contraindicates independence between D & E.")
            else:
                print("=> ChiSquared Test indicates independence between D & E.")
            print("")
        return pval
    def display(self):
        print("2x2 Contingency Table:")
        display(self.cont_table)
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
class Strata:
    def __init__(self, tables):
        self.tables = tables;
    def pool_table(self):
        a = sum([i.a for i in self.tables])
        b = sum([i.b for i in self.tables])
        c = sum([i.c for i in self.tables])
        d = sum([i.d for i in self.tables])
        return ContTable(a, b, c, d, 0)
    def cmh_test(self, verbose=1):
        Ai = [(i.a + i.b)*(i.a + i.c)/i.num_tot for i in self.tables]
        Vi = [(i.a + i.b)*(i.c + i.d)*(i.a + i.c)*(i.b + i.d)/(i.num_tot**2)/(i.num_tot-1) for i in self.tables]
    def woolf_correction(self, CI=.95, verbose=1):
        weights = [1/((1/(i.a+.5)) + (1/(i.b+.5)) + (1/(i.c+.5)) + (1/(i.d+.5)))  for i in self.tables]
        adj_logOR = sum([w*i.ss_odds_ratio(verbose=0)["logOR"] for i, w in zip(self.tables, weights)])/sum(weights)
        adj_OR = exp(adj_logOR)
        adj_logORVar = 1/sum(weights)
        z = -norm.interval(CI)[0]
        upper = exp(adj_logOR + z * sqrt(adj_logORVar))
        #Lower bound
        lower = exp(adj_logOR - z * sqrt(adj_logORVar))
        if verbose > 1:
            print("Individual ORs:", [i.odds_ratio(verbose=0)["OR"] for i in tables])
        if verbose:
            print("Woolf adj. OR:", adj_OR)
            if CI:
                print("=> {:d}% CI using logOR: {:f}".format(floor(CI*100), adj_logOR)) 
                print("=> => logOR Var:", adj_logORVar)
                print("=> => Upper Bound:", upper)
                print("=> => Lower Bound:", lower)
        return adj_OR, adj_logOR, (lower, upper)     
    def mantel_haenszel_correction(self, CI=.95, verbose=1):
        weights = [i.b*i.c/i.num_tot for i in self.tables]
        adj_OR = sum([i.a*i.d/i.num_tot for i in self.tables])/sum(weights)
        adj_logOR = log(adj_OR)

        first_num = sum([((i.a + i.d)/i.num_tot)*((i.a*i.d)/i.num_tot) for i in self.tables])
        first_denom = 2*sum([((i.a*i.d)/i.num_tot) for i in self.tables])**2
        second_num = sum([((i.a + i.d)/i.num_tot)*((i.b*i.c)/i.num_tot) + ((i.b + i.c)/i.num_tot)*((i.a*i.d)/i.num_tot) for i in self.tables])
        second_denom = 2*sum([((i.a*i.d)/i.num_tot) for i in self.tables])*sum([((i.b*i.c)/i.num_tot) for i in self.tables])
        third_num = sum([((i.b + i.c)/i.num_tot)*((i.b*i.c)/i.num_tot) for i in self.tables])
        third_denom = 2*sum([((i.b*i.c)/i.num_tot) for i in self.tables])**2
        adj_logORVar = (first_num/first_denom) + (second_num/second_denom) + (third_num/third_denom);

        z = -norm.interval(CI)[0]
        upper = exp(adj_logOR + z * sqrt(adj_logORVar))
        #Lower bound
        lower = exp(adj_logOR - z * sqrt(adj_logORVar))
        if verbose > 1:
            print("Individual ORs:", [i.odds_ratio(verbose=0)["OR"] for i in self.tables])
        if verbose:
            print("Mantel-Haenszel adj. OR:", adj_OR)
            if CI:
                print("=> {:d}% CI using logOR: {:f}".format(floor(CI*100), adj_logOR)) 
                print("=> => logOR Var:", adj_logORVar)
                print("=> => Upper Bound:", upper)
                print("=> => Lower Bound:", lower)
        return adj_OR, adj_logOR, (lower, upper)
