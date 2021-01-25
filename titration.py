import math
from sigfig import round
import decimal

soln_type = ''
titration_type = ''
eq_constant = ''
eq_type = ''
soln_conc = 0.0
soln_vol = 0.0
titrant_conc = 0.0
titrant_vol = 0.0
total_vol = 0.0
eq_vol = 0.0
vol_at_eq = 0.0
sig_figs = 0

# Main
def main():
    global soln_type, titration_type, eq_constant, eq_type, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, eq_vol, vol_at_eq, sig_figs

    soln_type = input("Type of Solution: ")
    titration_type = input("Type of Titration: ")    

    if (titration_type == 'weak-strong'):   
        eq_constant = input("Equilibrium Constant: ")
        eq_type = input("Ka or Kb? ")
        soln_conc  = input("Solution Concentration: ")
        soln_vol = input("Solution Volume: ")
        titrant_conc = input("Titrant Concentration: ")
        titrant_vol = input("Titrant Volume: ")
        sig_figs = input("Max Number of Sig Figs: ")
        eq_constant = float(eq_constant)
    elif (titration_type == 'strong-strong'):  
        soln_conc  = input("Solution Concentration: ")
        soln_vol = input("Solution Volume: ")
        titrant_conc = input("Titrant Concentration: ")
        titrant_vol = input("Titrant Volume: ")
        sig_figs = input("Max Number of Sig Figs: ")

    
    soln_conc = float(soln_conc)
    soln_vol = float(soln_vol)
    titrant_conc = float(titrant_conc)
    titrant_vol = float(titrant_vol)
    total_vol = float(total_vol)
    eq_vol = float(eq_vol)
    sig_figs = int(sig_figs)
    
    total_vol = soln_vol + titrant_vol
    eq_vol = (soln_conc * soln_vol) / titrant_conc
    eq_vol = round(eq_vol, sigfigs=sig_figs)
    vol_at_eq = eq_vol + soln_vol

    halfway_point = eq_vol / 2

    if (titration_type == 'weak-strong'):   
        if (titrant_vol < eq_vol):
            print("\n\nThe pH before the addition of any titrant is ", no_titrant_weak(eq_type, eq_constant, soln_conc, sig_figs))
            print("The pH before equivalence point is ", before_eqv_weak(eq_type, eq_constant, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs))
            print("The halfway point to Equivalence will be reached at ", halfway_point, " L of titrant where the pH is ", halfway_eqv_weak(eq_type, eq_constant, sig_figs))
            print("The equilvalence point will be reached at ", eq_vol, " L of titrant where the pH is ", at_eqv_weak(eq_type, eq_constant, titrant_conc, vol_at_eq, eq_vol, sig_figs), "\n\n\n")
        elif (titrant_vol > eq_vol):
            print("\n\nThe pH before the addition of any titrant is ", no_titrant_weak(eq_type, eq_constant, soln_conc, sig_figs))
            print("The halfway point to Equivalence will be reached at ", halfway_point, " L of titrant where the pH is ", halfway_eqv_weak(eq_type, eq_constant, sig_figs))
            print("The equilvalence point will be reached at ", eq_vol, " L of titrant where the pH is ", at_eqv_weak(eq_type, eq_constant, titrant_conc, vol_at_eq, eq_vol, sig_figs))
            print("The pH after equivalence point is ", after_eqv_weak(eq_type, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs), "\n\n\n")
        elif (titrant_vol == eq_vol):
            print("\n\n\nThis is the equivalence point. The pH is ", at_eqv_weak(eq_type, eq_constant, titrant_conc, vol_at_eq, eq_vol, sig_figs), "\n\n\n")
    elif (titration_type == 'strong-strong'):
        if (titrant_vol < eq_vol):
            print("\n\n\nThe pH before the addition of any titrant is", no_titrant_strong(soln_type, soln_conc, sig_figs))
            print("The pH before the equivlanece point is ", before_eqv_strong(soln_type, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs))
            print("The pH at the equivalence point is ", at_eqv_strong(), "\n\n\n")
        elif (titrant_vol > eq_vol):
            print("\n\n\nThe pH before the addition of any titrant is", no_titrant_strong(soln_type, soln_conc, sig_figs))
            print("The pH at the equivalence point is ", at_eqv_strong())
            print("The pH after the equivalence point is ", after_eqv_strong(soln_type, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs), "\n\n\n")
        elif (titrant_vol == eq_vol):
            print("\n\n\nThis is the equivalence point. The pH is ", at_eqv_strong(), "\n\n\n")

# Calculates the -log of numbers
def p_thing(num, sig_figs):
    ans = round((math.log10(num)) * (-1), decimals=sig_figs)
    return ans

# pOH to pH Conversion
def pOH_to_pH(concentration, sig_figs):
    pOH = p_thing(concentration, sig_figs)
    pH = 14.000000 - pOH
    pH = round(pH, decimals=sig_figs)
    return pH

# Swaps the Type of Equilibrium Constant
def eq_const_swap(eq_constant, sig_figs):
    converted_constant = (1e-14)/float(eq_constant)
    converted_constant = round(converted_constant, sigfigs=sig_figs)
    return converted_constant

# Henderson-Hasselbalch Equation
def henderson_hasselbalch(eq_type, eq_constant, base_c, acid_c, sig_figs):
    if (eq_type != 'Ka'):
        eq_const_swap(eq_constant, sig_figs)
    ratio = base_c / acid_c
    ans = p_thing(eq_constant, sig_figs) + ((p_thing(ratio, sig_figs)) * (-1))
    ans = round(ans, decimals=sig_figs)
    return ans

# Weak-Strong Titration: No Titrant Added         
def no_titrant_weak(eq_type, eq_constant, soln_conc, sig_figs):
    concentration = math.sqrt(eq_constant * soln_conc)
    concentration = round(concentration, sigfigs=sig_figs)
    if (eq_type == 'Ka'):
        ans = p_thing(concentration, sig_figs)
        return ans
    else:
        return pOH_to_pH(concentration, sig_figs)

# Weak-Strong Titration: Before Equivalence Point
def before_eqv_weak(eq_type, eq_constant, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs):
    titrant_mol = titrant_vol * titrant_conc
    soln_mol = (soln_vol * soln_conc) - titrant_mol
    titrant_M = titrant_mol / total_vol
    soln_M = soln_mol / total_vol
    if (eq_type == 'Ka'):
        return henderson_hasselbalch(eq_type, eq_constant, titrant_M, soln_M, sig_figs)
    else:
        cvt_constant = eq_const_swap(eq_constant, sig_figs)
        return henderson_hasselbalch(eq_type, cvt_constant, soln_M, titrant_M, sig_figs)
    
# Weak-Strong Titration: Halfway to Equivalence    
def halfway_eqv_weak(eq_type, eq_constant, sig_figs):
    if (eq_type == 'Ka'):
        return p_thing(eq_constant, sig_figs)
    else:
        cvt_constant = eq_const_swap(eq_constant, sig_figs)
        pH = p_thing(cvt_constant, sig_figs)
        return pH

# Weak-Strong Titration: At Equivalence
def at_eqv_weak(eq_type, eq_constant, titrant_conc, vol_at_eq, eq_vol, sig_figs):
    cvt_constant = eq_const_swap(eq_constant, sig_figs)
    rice_M = (eq_vol * titrant_conc) / vol_at_eq
    concentration = math.sqrt(cvt_constant * rice_M)
    concentration = round(concentration, sigfigs=sig_figs)
    if (eq_type == 'Ka'):
        return pOH_to_pH(concentration, sig_figs)
    else: 
        pH = p_thing(concentration, sig_figs)
        return pH #CHANGE TO pH

# General Use: After Equivalence Point
def after_eqv_weak(eq_type, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs):
    soln_mol = soln_vol * soln_conc
    titrant_mol = (titrant_vol * titrant_conc) - soln_mol
    titrant_M = titrant_mol / total_vol
    if (eq_type == 'Ka'):  
        return pOH_to_pH(titrant_M, sig_figs)
    else:
        return p_thing(titrant_M, sig_figs)
          
# Strong-Strong Titration: No Titrant Added
def no_titrant_strong(titration_type, soln_conc, sig_figs):
    if (titration_type == 'acidic'):
        return p_thing(soln_conc, sig_figs)
    else:
        return pOH_to_pH(soln_conc, sig_figs)

# Strong-Strong Titration: Before Equivalence Point
def before_eqv_strong(titration_type, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs):
    titrant_mol = titrant_vol * titrant_conc
    soln_mol = (soln_vol * soln_conc) - titrant_mol
    soln_M = soln_mol / total_vol
    if (titration_type == 'acidic'):
        return p_thing(soln_M, sig_figs)
    else:
        return pOH_to_pH(soln_M, sig_figs)

# Strong-Strong Titration: At Equivalence Point
def at_eqv_strong():
    return 7.00

# Strong-Strong Titration: After Equivalence Point
def after_eqv_strong(soln_type, soln_conc, soln_vol, titrant_conc, titrant_vol, total_vol, sig_figs):
    soln_mol = soln_vol * soln_conc
    titrant_mol = (titrant_vol * titrant_conc) - soln_mol
    titrant_M = titrant_mol / total_vol
    if (soln_type == 'acidic'):  
        return pOH_to_pH(titrant_M, sig_figs)
    else:
        return p_thing(titrant_M, sig_figs)

# MAIN BLOCK DO NOT TOUCH
if __name__ == '__main__':
    main()