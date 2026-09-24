import math
import itertools
try:
    from tabulate import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False

# ============================================================
#  1. Element Database  (r=pm, B=GPa, alpha=1/K, Tm=K)
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
#  2. Physics & DFT helpers
# ============================================================
def r_bar(els, fracs): return sum(fracs[i]*ELEMENTS_DATA[els[i]]["r"] for i in range(len(els)))
def delta_r(els, fracs):
    rb = r_bar(els, fracs)
    return math.sqrt(sum(fracs[i]*(1-ELEMENTS_DATA[els[i]]["r"]/rb)**2 for i in range(len(els))))*100
def delta_S_mix(fracs): return -R*sum(c*math.log(c) for c in fracs if c > 0)
def delta_H_mix(els, fracs):
    n = len(els)
    total = 0
    for i in range(n):
        for j in range(n):
            if i != j: total += 4*DELTA_H_IJ.get((els[i],els[j]),0)*fracs[i]*fracs[j]
    return total/2
def VEC_mix(els, fracs): return sum(fracs[i]*ELEMENTS_DATA[els[i]]["VEC"] for i in range(len(els)))
def Tm_mean(els, fracs): return sum(fracs[i]*ELEMENTS_DATA[els[i]]["Tm"] for i in range(len(els)))
def Omega(Tm, dS, dH): return float('inf') if abs(dH) < 1e-10 else Tm*dS/(abs(dH)*1000)
def gamma(els, fracs):
    rb = r_bar(els, fracs)
    rv = [ELEMENTS_DATA[el]["r"] for el in els]
    rmin, rmax = min(rv), max(rv)
    den = 1-(rmax/(rmax+rb))**2
    return float('inf') if den == 0 else (1-(rmin/(rmin+rb))**2)/den
def delta_chi(els, fracs):
    chi_bar = sum(fracs[i]*ELEMENTS_DATA[els[i]]["chi"] for i in range(len(els)))
    return math.sqrt(sum(fracs[i]*(1-ELEMENTS_DATA[els[i]]["chi"]/chi_bar)**2 for i in range(len(els))))*100
def ea_mix(els, fracs): return sum(fracs[i]*ELEMENTS_DATA[els[i]]["ea"] for i in range(len(els)))

# DFT Energy Mix (Rule of Mixtures for base energies)
def E_mix(els, fracs, phase):
    return sum(fracs[i]*ELEMENTS_DATA[els[i]][f"E_{phase}"] for i in range(len(els)))

# ============================================================
#  3. Criteria
# ============================================================
def check_SS(dH, dr, Om, gm, dchi, VEC, ea, T, Tm_val):
    results = {}
    results["MC1 Zhang (2008)"]   = "✅ SS" if (0.5<dr<6.5) and (-17.5<dH<5)     else "❌ IM/BMG"
    results["MC2 Yang (2012)"]    = "✅ SS" if (Om>=1.1) and (dr<=6.6)            else "❌ IM"
    results["MC3 Guo (2013)"]     = "✅ SS" if (-11.6<dH<3.2) and (dr<=6.6)      else "❌ IM/BMG"
    results["MC4 Wang γ (2014)"]  = "✅ SS" if gm <= 1.175                        else "❌ --"
    results["MC5 Poletti (2014)"] = "✅ SS" if (1<dr<6) and (3<dchi<6)            else "❌ --"
    ratio = T / Tm_val
    if ratio > 0.9:
        results["MC7 Wang (2014)"] = "✅ SS" if (-15<=dH<=5) and (dr<=6.6) else "❌ --"
    elif 0.5 <= ratio < 0.9:
        results["MC7 Wang (2014)"] = "✅ SS" if (dH>=-7.5) and (dr<=3.3) else "❌ --"
    else:
        results["MC7 Wang (2014)"] = "⚠️ T/Tm out of range"
    return results

def check_FCC(VEC, ea, T, Tm_val):
    results = {}
    # LSS1 — Guo (2011)
    if VEC >= 8: results["LSS1 Guo (2011)"] = "✅ FCC  (VEC≥8)"
    elif VEC < 6.87: results["LSS1 Guo (2011)"] = "❌ BCC  (VEC<6.87)"
    else: results["LSS1 Guo (2011)"] = "⚠️ Mixed"
    
    # LSS2 — Poletti (2014)
    if VEC > 7.5 and (1.6 < ea < 1.8): results["LSS2 Poletti (2014)"] = "✅ FCC"
    elif VEC < 7.5 and (1.8 < ea < 2.3): results["LSS2 Poletti (2014)"] = "❌ BCC"
    else: results["LSS2 Poletti (2014)"] = "⚠️ Unresolved"

    # LSS3 — Wang (2014)
    ratio = T / Tm_val
    if ratio > 0.9:
        results["LSS3 Wang (2014)"] = "✅ FCC" if VEC > 7.84 else ("❌ BCC" if VEC < 6.87 else "⚠️ Mixed")
    elif 0.5 < ratio < 0.9:
        results["LSS3 Wang (2014)"] = "✅ FCC" if VEC > 7.8 else ("❌ BCC" if VEC < 6 else "⚠️ Mixed")
    else:
        results["LSS3 Wang (2014)"] = "⚠️ T/Tm out of range"
        
    return results

# ============================================================
#  4. Generators & Formatters
# ============================================================
def generate_compositions(step, num_elements):
    steps = int(round(1.0 / step))
    for c in itertools.combinations_with_replacement(range(steps + 1), num_elements):
        if sum(c) == steps:
            for p in set(itertools.permutations(c)):
                if all(val > 0 for val in p):
                    yield tuple(val * step for val in p)

def show_criteria(title, criteria):
    print(f"\n  {'─'*6} {title} {'─'*(58-len(title))}")
    for key, val in criteria.items():
        marker = "  ✅" if "✅" in str(val) else ("  ❌" if "❌" in str(val) else "  ⚠️")
        print(f"{marker}  {key:<34} {val}")

def calculate_all_params(elements, fracs, T):
    dH   = delta_H_mix(elements, fracs)
    dS   = delta_S_mix(fracs)
    dr   = delta_r(elements, fracs)
    VEC  = VEC_mix(elements, fracs)
    gm   = gamma(elements, fracs)
    dchi = delta_chi(elements, fracs)
    ea   = ea_mix(elements, fracs)
    Tm   = Tm_mean(elements, fracs)
    Om   = Omega(Tm, dS, dH)
    e_fcc = E_mix(elements, fracs, "FCC")
    e_bcc = E_mix(elements, fracs, "BCC")
    e_hcp = E_mix(elements, fracs, "HCP")
    
    ss_results = check_SS(dH, dr, Om, gm, dchi, VEC, ea, T, Tm)
    fcc_results = check_FCC(VEC, ea, T, Tm)
    ss_pass = sum(1 for v in ss_results.values() if "✅" in v)
    fcc_pass = sum(1 for v in fcc_results.values() if "✅" in v)
    
    return {
        "fracs": fracs, "dH": dH, "dS": dS, "dr": dr, "VEC": VEC, "gm": gm, "dchi": dchi, 
        "ea": ea, "Tm": Tm, "Om": Om, "E_FCC": e_fcc, "E_BCC": e_bcc, "E_HCP": e_hcp,
        "ss_results": ss_results, "fcc_results": fcc_results,
        "ss_pass": ss_pass, "fcc_pass": fcc_pass
    }

# ============================================================
#  5. Main Logic
# ============================================================
def main():
    print("=" * 68)
    print("   FCC + SS RANGE CALCULATOR & DFT ENERGIES")
    print("=" * 68)
    
    raw = input("\n  Enter elements separated by space (e.g. Fe Co Ni Cr Mn): ")
    elements = [e.capitalize() for e in raw.split()]
    
    bad_els = [e for e in elements if e not in ELEMENTS_DATA]
    if bad_els:
        print(f"  ❌ Unknown elements: {bad_els}")
        return
        
    try: T = float(input("  Temperature (K) [e.g. 1000]: "))
    except ValueError: T = 1000.0

    step = 0.05
    min_ss_score = 3
    min_fcc_score = 1
    
    valid_data = []
    print(f"\n  ⚙️ Scanning combinations... (Step={step*100:.1f}%)")
    
    for fracs in generate_compositions(step, len(elements)):
        data = calculate_all_params(elements, fracs, T)
        # اعمال فیلتر (شروط قابل قبول)
        if data["ss_pass"] >= min_ss_score and data["fcc_pass"] >= min_fcc_score:
            valid_data.append(data)
            
    if not valid_data:
        print("\n  ❌ No combinations met the SS+FCC criteria.")
        return

    print(f"  ✅ Found {len(valid_data)} valid combinations.")
    
    # ── 1. Concentration Ranges ──
    print(f"\n{'='*68}\n  CONCENTRATION RANGES\n{'='*68}")
    print(f"  (Based on passing ≥{min_ss_score} SS criteria and ≥{min_fcc_score} FCC criteria)")
    for i, el in enumerate(elements):
        el_fracs = [d["fracs"][i]*100 for d in valid_data]
        print(f"    {el:2} :  {min(el_fracs):4.1f}%  to  {max(el_fracs):4.1f}%")

    # ── 2. Parameter Ranges ──
    print(f"\n{'='*68}\n  THERMODYNAMIC & DFT RANGES (MIN - MAX)\n{'='*68}")
    params_to_show = ["dH", "dr", "VEC", "Om", "gm", "dchi", "ea", "E_FCC", "E_BCC"]
    labels = {
        "dH": "ΔH_mix (kJ/mol)", "dr": "δ_r (%)", "VEC": "VEC (-)", 
        "Om": "Ω (Omega)", "gm": "γ (gamma)", "dchi": "Δχ (%)", 
        "ea": "e/a (-)", "E_FCC": "DFT E_FCC (eV)", "E_BCC": "DFT E_BCC (eV)"
    }
    
    for p in params_to_show:
        vals = [d[p] for d in valid_data]
        min_v, max_v = min(vals), max(vals)
        print(f"    {labels[p]:<20}:  {min_v:7.3f}   to   {max_v:7.3f}")

    # ── 3. The "Best" Example ──
    # پایدارترین FCC را بر اساس پایین‌ترین (منفی‌ترین) انرژی E_FCC یا بالاترین VEC انتخاب می‌کنیم
    best_data = sorted(valid_data, key=lambda x: x["E_FCC"])[0]  
    
    print(f"\n{'='*68}\n  EXAMPLE: MOST STABLE FCC COMPOSITION (Lowest E_FCC)\n{'='*68}")
    comp_str = " ".join([f"{el}{int(best_data['fracs'][i]*100)}" for i, el in enumerate(elements)])
    print(f"  Composition : {comp_str}")
    
    # چاپ جدول پارامترهای فرمول ها
    params_table = {
        "ΔH_mix  (kJ/mol)" : best_data["dH"],
        "ΔS_mix  (J/mol·K)": best_data["dS"],
        "δ_r  (%)"         : best_data["dr"],
        "VEC"              : best_data["VEC"],
        "Ω (Omega)"        : best_data["Om"],
        "γ (gamma)"        : best_data["gm"],
        "Δχ  (%)"          : best_data["dchi"],
        "e/a"              : best_data["ea"],
        "DFT E_FCC (eV)"   : best_data["E_FCC"],
        "DFT E_BCC (eV)"   : best_data["E_BCC"],
        "DFT E_HCP (eV)"   : best_data["E_HCP"],
    }
    
    print("\n  Calculated Values:")
    if HAS_TABULATE:
        rows = [(k, f"{v:.4f}") for k, v in params_table.items()]
        print(tabulate(rows, headers=["Parameter", "Value"], tablefmt="fancy_grid"))
    else:
        for k, v in params_table.items():
            print(f"    {k:<20}: {v:.4f}")
            
    show_criteria("Passed SS Criteria for this example", best_data["ss_results"])
    show_criteria("Passed FCC Criteria for this example", best_data["fcc_results"])
    print("=" * 68)

if __name__ == "__main__":
    main()