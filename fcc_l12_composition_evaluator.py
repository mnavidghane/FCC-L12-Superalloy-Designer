import math
try:
    from tabulate import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False

# ============================================================
#  1. Element Database & Enthalpy Matrix
# ============================================================
ELEMENTS_DATA = {
    "Nb": {"r":146,"VEC":5, "ea":1.51,"chi":1.60,"Tm":2477+273.15,"E_BCC":-10.0466,"E_FCC":-9.7232, "E_HCP":-9.7551},
    "Mo": {"r":139,"VEC":6, "ea":1.87,"chi":2.16,"Tm":2623+273.15,"E_BCC":-10.7799,"E_FCC":-10.3784,"E_HCP":-10.3666},
    "Ta": {"r":146,"VEC":5, "ea":1.50,"chi":1.50,"Tm":3017+273.15,"E_BCC":-11.7358,"E_FCC":-11.1896,"E_HCP":-11.4579},
    "W":  {"r":139,"VEC":6, "ea":1.86,"chi":2.36,"Tm":3422+273.15,"E_BCC":-12.7781,"E_FCC":-12.3115,"E_HCP":-12.2928},
    "V":  {"r":135,"VEC":5, "ea":1.50,"chi":1.63,"Tm":1910+273.15,"E_BCC":-8.9632, "E_FCC":-8.7150, "E_HCP":-8.7095},
    "Cr": {"r":128,"VEC":6, "ea":1.87,"chi":1.66,"Tm":1907+273.15,"E_BCC":-9.4655, "E_FCC":-9.0845, "E_HCP":-9.0751},
    "Ti": {"r":147,"VEC":4, "ea":1.33,"chi":1.54,"Tm":1668+273.15,"E_BCC":-2.2301, "E_FCC":-2.2155, "E_HCP":-2.2343},
    "Zr": {"r":160,"VEC":4, "ea":1.33,"chi":1.33,"Tm":1855+273.15,"E_BCC":-8.3598, "E_FCC":-8.3972, "E_HCP":-8.4354},
    "Hf": {"r":159,"VEC":4, "ea":1.32,"chi":1.30,"Tm":2233+273.15,"E_BCC":-9.6562, "E_FCC":-9.7613, "E_HCP":-9.8320},
    "Ni": {"r":125,"VEC":10,"ea":1.75,"chi":1.91,"Tm":1455+273.15,"E_BCC":-5.2954, "E_FCC":-5.3902, "E_HCP":-5.3681},
    "Co": {"r":125,"VEC":9, "ea":1.72,"chi":1.88,"Tm":1495+273.15,"E_BCC":-6.8834, "E_FCC":-6.9696, "E_HCP":-6.9902},
    "Fe": {"r":126,"VEC":8, "ea":1.83,"chi":1.83,"Tm":1538+273.15,"E_BCC":-8.2748, "E_FCC":-8.1872, "E_HCP":-8.2676},
    "Cu": {"r":128,"VEC":11,"ea":1.85,"chi":1.90,"Tm":1085+273.15,"E_BCC":-3.6082, "E_FCC":-3.6372, "E_HCP":-3.6325},
    "Mn": {"r":127,"VEC":7, "ea":1.61,"chi":1.55,"Tm":1246+273.15,"E_BCC":-8.8117, "E_FCC":-8.8885, "E_HCP":-8.9197},
    "Al": {"r":143,"VEC":3, "ea":1.36,"chi":1.61,"Tm": 660+273.15,"E_BCC":-3.6012, "E_FCC":-3.6967, "E_HCP":-3.6672},
}

DELTA_H_IJ = {
    # ── Refractory pairs ──
    ("Nb","Mo"):-6, ("Mo","Nb"):-6, ("Nb","Ta"): 0, ("Ta","Nb"): 0, ("Nb","W") :-8, ("W", "Nb"):-8,
    ("Nb","V") :-1, ("V", "Nb"):-1, ("Nb","Cr"):-7, ("Cr","Nb"):-7, ("Nb","Ti"): 2, ("Ti","Nb"): 2,
    ("Nb","Zr"): 4, ("Zr","Nb"): 4, ("Nb","Hf"): 4, ("Hf","Nb"): 4, ("Mo","Ta"):-5, ("Ta","Mo"):-5,
    ("Mo","W") : 0, ("W", "Mo"): 0, ("Mo","V") : 0, ("V", "Mo"): 0, ("Mo","Cr"): 0, ("Cr","Mo"): 0,
    ("Mo","Ti"):-4, ("Ti","Mo"):-4, ("Mo","Zr"):-6, ("Zr","Mo"):-6, ("Mo","Hf"):-4, ("Hf","Mo"):-4,
    ("Ta","W") :-7, ("W", "Ta"):-7, ("Ta","V") :-1, ("V", "Ta"):-1, ("Ta","Cr"):-7, ("Cr","Ta"):-7,
    ("Ta","Ti"): 1, ("Ti","Ta"): 1, ("Ta","Zr"): 3, ("Zr","Ta"): 3, ("Ta","Hf"): 3, ("Hf","Ta"): 3,
    ("W", "V") :-1, ("V", "W") :-1, ("W", "Cr"): 1, ("Cr","W") : 1, ("W", "Ti"):-6, ("Ti","W") :-6,
    ("W", "Zr"):-9, ("Zr","W") :-9, ("W", "Hf"):-6, ("Hf","W") :-6, ("V", "Cr"):-2, ("Cr","V") :-2,
    ("V", "Ti"):-2, ("Ti","V") :-2, ("V", "Zr"):-4, ("Zr","V") :-4, ("V", "Hf"):-2, ("Hf","V") :-2,
    ("Cr","Ti"):-7, ("Ti","Cr"):-7, ("Cr","Zr"):-12,("Zr","Cr"):-12,("Cr","Hf"):-9, ("Hf","Cr"):-9,
    ("Ti","Zr"): 0, ("Zr","Ti"): 0, ("Ti","Hf"): 0, ("Hf","Ti"): 0, ("Zr","Hf"): 0, ("Hf","Zr"): 0,
    # ── FCC element pairs ──
    ("Ni","Co"): 0, ("Co","Ni"): 0, ("Ni","Fe"):-2, ("Fe","Ni"):-2, ("Ni","Cu"): 4, ("Cu","Ni"): 4,
    ("Ni","Mn"):-8, ("Mn","Ni"):-8, ("Ni","Cr"):-7, ("Cr","Ni"):-7, ("Co","Fe"):-1, ("Fe","Co"):-1,
    ("Co","Cu"): 6, ("Cu","Co"): 6, ("Co","Mn"):-5, ("Mn","Co"):-5, ("Co","Cr"):-4, ("Cr","Co"):-4,
    ("Fe","Cu"):13, ("Cu","Fe"):13, ("Fe","Mn"): 0, ("Mn","Fe"): 0, ("Fe","Cr"):-1, ("Cr","Fe"):-1,
    ("Cu","Mn"): 4, ("Mn","Cu"): 4, ("Cu","Cr"):12, ("Cr","Cu"):12, ("Mn","Cr"): 2, ("Cr","Mn"): 2,
    # ── Al pairs ──
    ("Al","Ni"):-22,("Ni","Al"):-22,("Al","Co"):-19,("Co","Al"):-19,("Al","Fe"):-11,("Fe","Al"):-11,
    ("Al","Cu"):-1, ("Cu","Al"):-1, ("Al","Mn"):-19,("Mn","Al"):-19,("Al","Cr"):-10,("Cr","Al"):-10,
    ("Al","Nb"):-18,("Nb","Al"):-18,("Al","Mo"):-19,("Mo","Al"):-19,("Al","Ta"):-19,("Ta","Al"):-19,
    ("Al","W") :-9, ("W", "Al"):-9, ("Al","V") :-16,("V", "Al"):-16,("Al","Ti"):-30,("Ti","Al"):-30,
    ("Al","Zr"):-44,("Zr","Al"):-44,("Al","Hf"):-39,("Hf","Al"):-39,
    # ── FCC ↔ Refractory ──
    ("Ni","Nb"):-30,("Nb","Ni"):-30,("Ni","Mo"):-7, ("Mo","Ni"):-7, ("Ni","Ta"):-29,("Ta","Ni"):-29,
    ("Ni","W") :-3, ("W", "Ni"):-3, ("Ni","V") :-18,("V", "Ni"):-18,("Ni","Ti"):-35,("Ti","Ni"):-35,
    ("Ni","Zr"):-49,("Zr","Ni"):-49,("Ni","Hf"):-42,("Hf","Ni"):-42,("Co","Nb"):-25,("Nb","Co"):-25,
    ("Co","Mo"):-5, ("Mo","Co"):-5, ("Co","Ta"):-24,("Ta","Co"):-24,("Co","W") :-1, ("W", "Co"):-1,
    ("Co","V") :-14,("V", "Co"):-14,("Co","Ti"):-28,("Ti","Co"):-28,("Co","Zr"):-41,("Zr","Co"):-41,
    ("Co","Hf"):-37,("Hf","Co"):-37,("Fe","Nb"):-16,("Nb","Fe"):-16,("Fe","Mo"):-2, ("Mo","Fe"):-2,
    ("Fe","Ta"):-15,("Ta","Fe"):-15,("Fe","W") : 0, ("W", "Fe"): 0, ("Fe","V") :-7, ("V", "Fe"):-7,
    ("Fe","Ti"):-17,("Ti","Fe"):-17,("Fe","Zr"):-25,("Zr","Fe"):-25,("Fe","Hf"):-22,("Hf","Fe"):-22,
    ("Cu","Nb"): 3, ("Nb","Cu"): 3, ("Cu","Mo"):19, ("Mo","Cu"):19, ("Cu","Ta"): 2, ("Ta","Cu"): 2,
    ("Cu","W") :22, ("W", "Cu"):22, ("Cu","V") : 5, ("V", "Cu"): 5, ("Cu","Ti"):-9, ("Ti","Cu"):-9,
    ("Cu","Zr"):-23,("Zr","Cu"):-23,("Cu","Hf"):-17,("Hf","Cu"):-17,("Mn","Nb"):-4, ("Nb","Mn"):-4,
    ("Mn","Mo"):-5, ("Mo","Mn"):-5, ("Mn","Ta"):-7, ("Ta","Mn"):-7, ("Mn","W") :-4, ("W", "Mn"):-4,
    ("Mn","V") :-5, ("V", "Mn"):-5, ("Mn","Ti"):-8, ("Ti","Mn"):-8, ("Mn","Zr"):-10,("Zr","Mn"):-10,
    ("Mn","Hf"):-8, ("Hf","Mn"):-8,
}
R = 8.314

# ============================================================
#  2. Physics Helpers
# ============================================================
def delta_S_mix(fracs): 
    return -R*sum(c*math.log(c) for c in fracs if c > 0)

def delta_H_mix(els, fracs):
    total = 0
    n = len(els)
    for i in range(n):
        for j in range(n):
            if i != j: total += 4*DELTA_H_IJ.get((els[i],els[j]),0)*fracs[i]*fracs[j]
    return total/2

def VEC_mix(els, fracs): 
    return sum(fracs[i]*ELEMENTS_DATA[els[i]]["VEC"] for i in range(len(els)))

def Tm_mean(els, fracs): 
    return sum(fracs[i]*ELEMENTS_DATA[els[i]]["Tm"] for i in range(len(els)))

def delta_chi_absolute(els, fracs):
    chi_bar = sum(fracs[i]*ELEMENTS_DATA[els[i]]["chi"] for i in range(len(els)))
    return math.sqrt(sum(fracs[i]*(ELEMENTS_DATA[els[i]]["chi"] - chi_bar)**2 for i in range(len(els))))

def get_element_percentage(els, fracs, target_el):
    if target_el in els:
        idx = els.index(target_el)
        return fracs[idx] * 100
    return 0.0

# ============================================================
#  3. Criteria Evaluator
# ============================================================
def evaluate_rules(els, fracs):
    dH = delta_H_mix(els, fracs)
    dS = delta_S_mix(fracs)
    VEC = VEC_mix(els, fracs)
    Tm = Tm_mean(els, fracs)
    dChi = delta_chi_absolute(els, fracs)
    
    co_pct = get_element_percentage(els, fracs, "Co")
    fe_pct = get_element_percentage(els, fracs, "Fe")
    al_pct = get_element_percentage(els, fracs, "Al")
    
    results = {}
    
    # ── 3 Main Core Rules (L12 Formation) ──
    results["Core: VEC > 8"] = ("✅" if VEC > 8 else "❌", f"VEC = {VEC:.2f}")
    results["Core: -16.0 < ΔH_mix < -9.7"] = ("✅" if -16.0 < dH < -9.7 else "❌", f"ΔH = {dH:.2f} kJ/mol")
    results["Core: 1671 < Tm < 1822 K"] = ("✅" if 1671 < Tm < 1822 else "❌", f"Tm = {Tm:.0f} K")
    
    # ── Refined Rules (Detrimental Phase Avoidance) ──
    results["Refined: 5.3 <= ΔS_mix <= 13.4"] = ("✅" if 5.3 <= dS <= 13.4 else "❌", f"ΔS = {dS:.2f} J/(mol·K)")
    results["Refined: Δχ < 0.12"] = ("✅" if dChi < 0.12 else "❌", f"Δχ = {dChi:.3f}")
    
    # ── Elemental Constraints (Based on the Paper) ──
    results["Constraint: Co content (31% to 72%)"] = ("✅" if 31 <= co_pct <= 72 else "❌", f"Co = {co_pct:.1f}%")
    results["Constraint: Fe content (<= 10%)"] = ("✅" if fe_pct <= 10 else "❌", f"Fe = {fe_pct:.1f}%")
    results["Constraint: Al content (<= 9%)"] = ("✅" if al_pct <= 9 else "❌", f"Al = {al_pct:.1f}%")
    
    all_passed = all(val[0] == "✅" for val in results.values())
    
    return all_passed, results, {"dH":dH, "dS":dS, "VEC":VEC, "Tm":Tm, "dChi":dChi}

# ============================================================
#  4. User Interface
# ============================================================
def ask_elements():
    print("\n" + "="*68)
    print("   EVALUATE YOUR ALLOY: FCC + L12 (DETRIMENTAL PHASE FILTER)")
    print("="*68)
    
    while True:
        raw = input("\n  Enter elements separated by space (e.g. Ni Co Cr Al Ti): ").strip()
        elements = [e.capitalize() for e in raw.split()]
        
        if not elements:
            continue
            
        bad_els = [e for e in elements if e not in ELEMENTS_DATA]
        if bad_els:
            print(f"  ❌ Unknown elements: {', '.join(bad_els)}")
            continue
            
        if len(set(elements)) != len(elements):
            print("  ❌ Duplicate elements detected.")
            continue
            
        return elements

def ask_composition(elements):
    print("\n  Enter atomic percentages (%) for each element (must sum to 100):")
    while True:
        vals = []
        for el in elements:
            while True:
                try:
                    v = float(input(f"    {el} (%): "))
                    if v < 0:
                        print("    ❌ Must be ≥ 0.")
                        continue
                    vals.append(v)
                    break
                except ValueError:
                    print("    ❌ Please enter a valid number.")
                    
        total = sum(vals)
        if abs(total - 100) > 0.5:
            print(f"  ❌ Sum is {total:.2f}%. It must be exactly 100%. Please try again.")
            continue
            
        # Normalize to exactly 1.0 (fractions)
        return [v / total for v in vals], [v / total * 100 for v in vals]

def print_results(elements, pcts, all_passed, results, params):
    print("\n" + "="*68)
    print("   ALLOY EVALUATION REPORT")
    print("="*68)
    
    comp_str = " ".join([f"{el}{pct:.1f}" for el, pct in zip(elements, pcts)])
    print(f"\n  🎯 Composition: {comp_str}")
    
    print("\n  📊 Thermodynamic Parameters:")
    param_rows = [
        ("VEC", f"{params['VEC']:.2f}"),
        ("ΔH_mix (kJ/mol)", f"{params['dH']:.2f}"),
        ("Tm (K)", f"{params['Tm']:.0f}"),
        ("ΔS_mix (J/mol·K)", f"{params['dS']:.2f}"),
        ("Δχ", f"{params['dChi']:.3f}")
    ]
    if HAS_TABULATE:
        print(tabulate(param_rows, headers=["Parameter", "Value"], tablefmt="fancy_grid"))
    else:
        for k, v in param_rows:
            print(f"    {k:<20}: {v}")

    print("\n  🔎 FCC + L12 Rule Evaluation:")
    rule_rows = []
    for rule, (status, val) in results.items():
        rule_rows.append((status, rule, val))
        
    if HAS_TABULATE:
        print(tabulate(rule_rows, headers=["Status", "Rule/Criterion", "Calculated"], tablefmt="fancy_grid"))
    else:
        for status, rule, val in rule_rows:
            print(f"    {status}  {rule:<35} | {val}")
            
    print("\n" + "="*68)
    if all_passed:
        print("   ✅ VERDICT: EXCELLENT! ")
        print("   This alloy perfectly satisfies all ML rules for 'FCC + L12'")
        print("   and effectively avoids known detrimental phases.")
    else:
        print("   ❌ VERDICT: FAILS ONE OR MORE CRITERIA.")
        print("   This alloy is unlikely to form a pure FCC + L12 dual-phase")
        print("   without detrimental phases. Please adjust the composition.")
    print("="*68 + "\n")

def main():
    while True:
        elements = ask_elements()
        fracs, pcts = ask_composition(elements)
        
        all_passed, results, params = evaluate_rules(elements, fracs)
        print_results(elements, pcts, all_passed, results, params)
        
        again = input("  Would you like to evaluate another composition? (y/n): ").strip().lower()
        if again != 'y':
            print("  Goodbye! 👋")
            break

if __name__ == "__main__":
    main()