import subprocess
import re

HARTREE_TO_EV = 1. #27.2114

def run_xtb(input_xyz, output_file='dipro_full.out'):
    with open(output_file, 'w') as out:
        subprocess.run(['xtb', input_xyz, '--dipro', '--gfn', '1'], stdout=out)

def extract_section(lines, header):
    """Extract block of data following a specific orbital combination header."""
    section = []
    in_section = False
    for line in lines:
        if header in line:
            in_section = True
            continue
        if in_section:
            if line.strip() == "":
                break
            section.append(line.strip())
    return section

def parse_block(block):
    """Extract E1, E2, J, S, and reported J_eff from a block of DIPRO output."""
    data = {}
    for line in block:
        if "E_mon(orb)" in line:
            parts = re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)
            data["E1"] = float(parts[-2])
            data["E2"] = float(parts[-1])
        elif "J(AB):" in line:
            data["J"] = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)[0])
        elif "S(AB):" in line:
            data["S"] = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)[0])
        elif "J(AB,eff):" in line and "J_eff_reported" not in data:
            data["J_eff_reported"] = float(re.findall(r"[-+]?[0-9]*\.?[0-9]+", line)[0])
    return data

def recalculate_J_eff(data):
    S_ev = data["S"] * HARTREE_TO_EV
    E_avg = (data["E1"] + data["E2"]) / 2
    numerator = data["J"] * 1.928 - S_ev * E_avg
    denominator = 1 - S_ev ** 2
    J_eff_calc = numerator / denominator
    return J_eff_calc, S_ev, E_avg

def parse_all(output_file):
    with open(output_file, 'r') as f:
        lines = f.readlines()

    headers = {
        "HOMO–HOMO": "Fragment 1 HOMO   Fragment 2 HOMO",
        "HOMO–LUMO": "Fragment 1 HOMO   Fragment 2 LUMO",
        "LUMO–HOMO": "Fragment 1 LUMO   Fragment 2 HOMO",
        "LUMO–LUMO": "Fragment 1 LUMO   Fragment 2 LUMO",
    }

    results = {}

    for key, header in headers.items():
        block = extract_section(lines, header)
        if block:
            data = parse_block(block)
            J_eff_calc, S_ev, E_avg = recalculate_J_eff(data)
            results[key] = {
                "E1": data["E1"],
                "E2": data["E2"],
                "J": data["J"],
                "S": data["S"],
                "S_ev": S_ev,
                "E_avg": E_avg,
                "J_eff_reported": data.get("J_eff_reported", None),
                "J_eff_calc": J_eff_calc,
            }

    return results

def summarize(results):
    print("\n=== DIPRO High-Precision Summary ===")
    total_CT = 0
    for label, res in results.items():
        print(f"\n--- {label} ---")
        print(f"E1:                   {res['E1']:12.6f} eV")
        print(f"E2:                   {res['E2']:12.6f} eV")
        print(f"J(AB):                {res['J']:12.6f} eV")
        print(f"S(AB):                {res['S']:12.8f} Eh ({res['S_ev']:.8f} eV)")
        print(f"E_avg:                {res['E_avg']:12.6f} eV")
        if res["J_eff_reported"] is not None:
            print(f"Reported J_eff:       {res['J_eff_reported']:12.6f} eV")
        print(f"Recalculated J_eff:   {res['J_eff_calc']:12.10f} eV")
        
        # Collect CT terms
        if label in ["HOMO–LUMO", "LUMO–HOMO"]:
            total_CT += res["J_eff_calc"]

    # Summary values
    print("\n=== Transport Summary ===")
    if "HOMO–HOMO" in results:
        print(f"Hole transport J_eff:        {results['HOMO–HOMO']['J_eff_calc']:.10f} eV")
    if "LUMO–LUMO" in results:
        print(f"Electron transport J_eff:    {results['LUMO–LUMO']['J_eff_calc']:.10f} eV")
    print(f"Charge transfer total J_eff: {total_CT:.10f} eV")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python dipro_full_parser.py input.xyz")
        sys.exit(1)

    input_xyz = sys.argv[1]
    run_xtb(input_xyz)
    results = parse_all("dipro_full.out")
    summarize(results)
