import math
import itertools
try:
    from tabulate import tabulate
    HAS_TABULATE = True
except ImportError:
    HAS_TABULATE = False

# ============================================================
#  1. Element Database
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
    ("Al","Ni"):-22,("Ni","Al"):-22,("Al","Co"):-19,("Co","Al"):-19,("Al","Fe"):-11,("Fe","Al"):-11,
    ("Al","Cr"):-10,("Cr","Al"):-10,("Al","Ti"):-30,("Ti","Al"):-30,("Ni","Co"): 0, ("Co","Ni"): 0,
    ("Ni","Fe"):-2, ("Fe","Ni"):-2, ("Ni","Cr"):-7, ("Cr","Ni"):-7, ("Co","Fe"):-1, ("Fe","Co"):-1,
    ("Co","Cr"):-4, ("Cr","Co"):-4, ("Fe","Cr"):-1, ("Cr","Fe"):-1, ("Ni","Ti"):-35,("Ti","Ni"):-35,
    ("Co","Ti"):-28,("Ti","Co"):-28,("Fe","Ti"):-17,("Ti","Fe"):-17,("Cr","Ti"):-7, ("Ti","Cr"):-7,
    # مقادیر سایر جفت‌ها در صورت لزوم از فایل قبلی خوانده می‌شود (برای سادگی اینجا خلاصه شده)
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
    """ محاسبه انحراف معیار مطلق الکترونگاتیویته (بدون درصد) همانطور که در مقاله آمده است """
    chi_bar = sum(fracs[i]*ELEMENTS_DATA[els[i]]["chi"] for i in range(len(els)))
    return math.sqrt(sum(fracs[i]*(ELEMENTS_DATA[els[i]]["chi"] - chi_bar)**2 for i in range(len(els))))

def E_mix(els, fracs, phase):
    return sum(fracs[i]*ELEMENTS_DATA[els[i]][f"E_{phase}"] for i in range(len(els)))

# ============================================================
#  3. Advanced ML Criteria (FCC + L12 + NO Other Phases)
# ============================================================
def check_advanced_rules(els, fracs, VEC, dH, Tm, dS, dChi):
    results = {}
    
    # ── 3 Main Core Rules ──
    results["Core 1: VEC > 8"] = "✅" if VEC > 8 else f"❌ ({VEC:.2f})"
    results["Core 2: -16.0 < ΔH_mix < -9.7"] = "✅" if -16.0 < dH < -9.7 else f"❌ ({dH:.2f})"
    results["Core 3: 1671 < Tm < 1822 K"] = "✅" if 1671 < Tm < 1822 else f"❌ ({Tm:.0f} K)"
    
    # ── Refined Rules (Detrimental Phase Avoidance) ──
    results["Refined 1: 5.3 <= ΔS_mix <= 13.4"] = "✅" if 5.3 <= dS <= 13.4 else f"❌ ({dS:.2f})"
    results["Refined 2: Δχ < 0.12"] = "✅" if dChi < 0.12 else f"❌ ({dChi:.3f})"
    
    # ── Elemental Constraints ──
    for i, el in enumerate(els):
        pct = fracs[i] * 100
        if el == "Co":
            results["Elem: Co (31% to 72%)"] = "✅" if 31 <= pct <= 72 else f"❌ ({pct:.1f}%)"
        if el == "Fe":
            results["Elem: Fe (<= 10%)"] = "✅" if pct <= 10 else f"❌ ({pct:.1f}%)"
        if el == "Al":
            results["Elem: Al (<= 9%)"] = "✅" if pct <= 9 else f"❌ ({pct:.1f}%)"
            
    return results

# ============================================================
#  4. Generators & Core Logic
# ============================================================
def generate_compositions(step, num_elements):
    steps = int(round(1.0 / step))
    for c in itertools.combinations_with_replacement(range(steps + 1), num_elements):
        if sum(c) == steps:
            for p in set(itertools.permutations(c)):
                if all(val > 0 for val in p):
                    yield tuple(val * step for val in p)

def evaluate_composition(els, fracs):
    dH = delta_H_mix(els, fracs)
    dS = delta_S_mix(fracs)
    VEC = VEC_mix(els, fracs)
    Tm = Tm_mean(els, fracs)
    dChi = delta_chi_absolute(els, fracs)
    
    rules_dict = check_advanced_rules(els, fracs, VEC, dH, Tm, dS, dChi)
    passed = all("✅" in v for v in rules_dict.values())
    
    return passed, rules_dict, {"dH":dH, "dS":dS, "VEC":VEC, "Tm":Tm, "dChi":dChi}

def main():
    print("=" * 68)
    print("   ADVANCED ML DESIGNER: FCC + L12 (DETRIMENTAL PHASE FILTER)")
    print("=" * 68)
    
    # معمولاً برای L12 این سیستم‌ها بسیار معروف هستند: Ni-Co-Cr-Al-Ti یا Ni-Co-Fe-Cr-Al
    raw = input("\n  Enter elements (e.g. Ni Co Cr Al Ti): ")
    elements = [e.capitalize() for e in raw.split()]
    
    bad_els = [e for e in elements if e not in ELEMENTS_DATA]
    if bad_els:
        print(f"  ❌ Unknown elements: {bad_els}"); return

    step = 0.02 # دقت جستجو را بالاتر بردم (استپ 2 درصدی) چون شروط سختگیرانه‌تر شده‌اند
    valid_data = []
    total_checked = 0
    
    print(f"\n  ⚙️ Scanning with strict rules... (Step={step*100:.0f}%)")
    
    for fracs in generate_compositions(step, len(elements)):
        total_checked += 1
        passed, rules, params = evaluate_composition(elements, fracs)
        if passed:
            params["fracs"] = fracs
            params["E_FCC"] = E_mix(elements, fracs, "FCC")
            params["rules"] = rules
            valid_data.append(params)
            
    print(f"  📊 Checked {total_checked} combinations.")
    
    if not valid_data:
        print("\n  ❌ No combinations met ALL the strict criteria.")
        print("  💡 Tips based on the paper:")
        print("     - Ensure Co is present and can reach 31-72%.")
        print("     - Ensure Fe can be kept under 10%.")
        print("     - Ensure Al can be kept under 9%.")
        return

    print(f"  ✅ Found {len(valid_data)} perfectly optimal 'FCC + L12' combinations.")
    
    # ── 1. Concentration Ranges ──
    print(f"\n{'='*68}\n  SAFE CONCENTRATION RANGES\n{'='*68}")
    for i, el in enumerate(elements):
        el_fracs = [d["fracs"][i]*100 for d in valid_data]
        print(f"    {el:2} :  {min(el_fracs):4.1f}%  to  {max(el_fracs):4.1f}%")

    # ── 2. Best Example ──
    best_data = sorted(valid_data, key=lambda x: x["E_FCC"])[0]  
    
    print(f"\n{'='*68}\n  EXAMPLE: BEST FCC+L12 CANDIDATE\n{'='*68}")
    comp_str = " ".join([f"{el}{best_data['fracs'][i]*100:.1f}" for i, el in enumerate(elements)])
    print(f"  Composition : {comp_str}\n")
    
    if HAS_TABULATE:
        rows = [(k, v) for k, v in best_data['rules'].items()]
        print(tabulate(rows, headers=["Criterion", "Status"], tablefmt="fancy_grid"))
    else:
        for k, v in best_data['rules'].items():
            print(f"    {k:<35}: {v}")

    print("=" * 68)

if __name__ == "__main__":
    main()