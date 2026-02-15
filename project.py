import re
import tkinter as tk
from tkinter import messagebox
from sympy import Matrix, lcm

SUBSCRIPT_MAP = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

VALID_ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Fl",
    "Lv", "Ts", "Og"
}

ATOMIC_MASSES = {
    "H": 1.008, "He": 4.003, "Li": 6.941, "Be": 9.012, "B": 10.81, "C": 12.01, "N": 14.01,
    "O": 16.00, "F": 19.00, "Ne": 20.18, "Na": 22.99, "Mg": 24.31, "Al": 26.98, "Si": 28.09,
    "P": 30.97, "S": 32.07, "Cl": 35.45, "Ar": 39.95, "K": 39.10, "Ca": 40.08, "Sc": 44.96,
    "Ti": 47.87, "V": 50.94, "Cr": 52.00, "Mn": 54.94, "Fe": 55.85, "Co": 58.93, "Ni": 58.69,
    "Cu": 63.55, "Zn": 65.38, "Ga": 69.72, "Ge": 72.63, "As": 74.92, "Se": 78.97, "Br": 79.90,
    "Kr": 83.80, "Rb": 85.47, "Sr": 87.62, "Y": 88.91, "Zr": 91.22, "Nb": 92.91, "Mo": 95.95,
    "Tc": 98.00, "Ru": 101.1, "Rh": 102.9, "Pd": 106.4, "Ag": 107.9, "Cd": 112.4, "In": 114.8,
    "Sn": 118.7, "Sb": 121.8, "Te": 127.6, "I": 126.9, "Xe": 131.3, "Cs": 132.9, "Ba": 137.3,
    "La": 138.9, "Ce": 140.1, "Pr": 140.9, "Nd": 144.2, "Pm": 145.0, "Sm": 150.4, "Eu": 152.0,
    "Gd": 157.3, "Tb": 158.9, "Dy": 162.5, "Ho": 164.9, "Er": 167.3, "Tm": 168.9, "Yb": 173.1,
    "Lu": 175.0, "Hf": 178.5, "Ta": 180.9, "W": 183.8, "Re": 186.2, "Os": 190.2, "Ir": 192.2,
    "Pt": 195.1, "Au": 197.0, "Hg": 200.6, "Tl": 204.4, "Pb": 207.2, "Bi": 209.0, "Po": 209.0,
    "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0, "Th": 232.0, "Pa": 231.0,
    "U": 238.0, "Np": 237.0, "Pu": 244.0, "Am": 243.0, "Cm": 247.0, "Bk": 247.0, "Cf": 251.0,
    "Es": 252.0, "Fm": 257.0, "Md": 258.0, "No": 259.0, "Lr": 262.0, "Rf": 267.0, "Db": 268.0,
    "Sg": 271.0, "Bh": 272.0, "Hs": 270.0, "Mt": 276.0, "Ds": 281.0, "Rg": 280.0, "Cn": 285.0,
    "Fl": 289.0, "Lv": 293.0, "Ts": 294.0, "Og": 294.0
}

def format_compound(compound):
    return compound.translate(SUBSCRIPT_MAP)

def check_parentheses(formula):
    stack = []
    for char in formula:
        if char == "(":
            stack.append(char)
        elif char == ")":
            if not stack:
                return False
            stack.pop()
    return len(stack) == 0

def parse_equation(reactants, products):
    reactants_list = [c.strip() for c in reactants.split("+") if c.strip()]
    products_list = [c.strip() for c in products.split("+") if c.strip()]
    
    if not reactants_list or not products_list:
        raise ValueError("equation must have reactants and products")
    
    return reactants_list, products_list

def count_atoms(formula):
    if not check_parentheses(formula):
        raise ValueError(f"missing parentheses in {formula}")
    
    stack = [{}]
    i = 0
    
    while i < len(formula):
        
        if formula[i] == "(":
            stack.append({})
            i += 1
        
        elif formula[i] == ")":
            i += 1
            num = ""
            while i < len(formula) and formula[i].isdigit():
                num += formula[i]
                i += 1
            
            multiplier = int(num) if num else 1
            top = stack.pop()
            
            for element in top:
                top[element] *= multiplier
            
            for element, count in top.items():
                stack[-1][element] = stack[-1].get(element, 0) + count
        
        else:
            match = re.match(r"^([A-Z][a-z]?)(\d*)", formula[i:])
            
            if not match:
                raise ValueError(f"invalid formula structure in {formula}")
            element, num = match.groups()
            
            if element not in VALID_ELEMENTS:
                raise ValueError(f"invalid element symbol: {element}")
            
            count = int(num) if num else 1
            stack[-1][element] = stack[-1].get(element, 0) + count
            
            i += len(match.group(0))
    
    return stack[0]

def calculate_molar_mass(formula):
    atoms = count_atoms(formula)
    return sum(ATOMIC_MASSES[element] * count for element, count in atoms.items())

def balance_equation(reactants, products):
    reactants_list, products_list = parse_equation(reactants, products)
    compounds = reactants_list + products_list
    compound_atom_counts = [count_atoms(c) for c in compounds]
    elements = set()
    
    for i, compound in enumerate(compounds):
        elements.update(compound_atom_counts[i].keys())
    
    elements = list(elements)
    matrix = []
    
    for element in elements:
        row = []
        for i in range(len(reactants_list)):
            row.append(compound_atom_counts[i].get(element, 0))
    
        for i in range(len(reactants_list), len(compounds)):
            row.append(-compound_atom_counts[i].get(element, 0))
    
        matrix.append(row)
    
    m = Matrix(matrix)
    nullspace = m.nullspace()
    
    if not nullspace:
        raise ValueError("can't be balanced - no solution :(")
    
    coeffs = nullspace[0]
    lcm_val = lcm([term.q for term in coeffs])
    coeffs = coeffs * lcm_val
    coeffs = [int(c) for c in coeffs]
    
    if all(c <= 0 for c in coeffs):
        coeffs = [-c for c in coeffs]
    
    if any(c == 0 for c in coeffs):
        raise ValueError("can't be balanced - zero coefficient found :(")
    
    result = []
    
    for i, compound in enumerate(compounds):
        coefficient = coeffs[i]
        formatted = format_compound(compound)
        result.append(f"{coefficient if coefficient != 1 else ''}{formatted}")
    balanced_eq_str = " + ".join(result[:len(reactants_list)]) + " -> " + " + ".join(result[len(reactants_list):])
    
    return balanced_eq_str, coeffs, compounds

class SubscriptEntry(tk.Entry):

    def __init__(self, *args, **kwargs):
        tk.Entry.__init__(self, *args, **kwargs)
        self.var = tk.StringVar()
        self.config(textvariable=self.var)
        self.var.trace_add("write", self.update_subscript)

    def update_subscript(self, *args):
        value = self.var.get()
        new_value = value.translate(SUBSCRIPT_MAP)
        if value != new_value:
            self.var.set(new_value)

class MainMenuGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🧪 chemistry toolkit 🧪")
        self.root.geometry("500x300")
        self.root.resizable(False, False)
        
        tk.Label(root, text="🧪 chemistry toolkit 🧪", font=("Consolas", 18, "bold")).pack(pady=30)
        tk.Label(root, text="choose your tool:", font=("Consolas", 12)).pack(pady=10)
        btn_frame = tk.Frame(root)
        btn_frame.pack(pady=20)
        
        tk.Button(btn_frame, text="equation balancer", font=("Consolas", 12), bg="#4CAF50", fg="white", width=30, height=2, command=self.open_balancer).pack(pady=10)
        
        tk.Button(btn_frame, text="stoichiometry calculator", font=("Consolas", 12), bg="#FF9800", fg="white", width=30, height=2, command=self.open_stoichiometry).pack(pady=10)
    
    def open_balancer(self):
        ChemBalancerGUI(tk.Toplevel(self.root))
    
    def open_stoichiometry(self):
        StoichiometryGUI(tk.Toplevel(self.root))

class ChemBalancerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🧪 equation balancer")
        self.root.geometry("650x250")
        self.root.resizable(False, False)
        
        tk.Label(root, text="enter your unbalanced chemical equation:", font=("Consolas", 12)).pack(pady=10)
        input_frame = tk.Frame(root)
        input_frame.pack(pady=5)
        
        react_frame = tk.Frame(input_frame)
        react_frame.grid(row=0, column=0, padx=10)
        tk.Label(react_frame, text="reactants", font=("Consolas", 10)).pack()
        self.react_entry = SubscriptEntry(react_frame, font=("Consolas", 14), width=25)
        
        self.react_entry.pack()
        
        tk.Label(input_frame, text="→", font=("Consolas", 20)).grid(row=0, column=1, padx=5)
        
        prod_frame = tk.Frame(input_frame)
        prod_frame.grid(row=0, column=2, padx=10)
        tk.Label(prod_frame, text="products", font=("Consolas", 10)).pack()
        self.prod_entry = SubscriptEntry(prod_frame, font=("Consolas", 14), width=25)
        self.prod_entry.pack()
        
        tk.Button(root, text="balance", font=("Consolas", 12), bg="#4CAF50", fg="white", command=self.balance).pack(pady=10)
        
        self.output = tk.Label(root, text="", font=("Consolas", 14), fg="blue", wraplength=600)
        self.output.pack(pady=20)
    
    def balance(self):
        reactants = self.react_entry.get().strip()
        products = self.prod_entry.get().strip()
        if not reactants or not products:
            messagebox.showwarning("hey!", "you gotta enter both reactants and products!")
            return
        try:
            balanced, _, _ = balance_equation(reactants, products)
            self.output.config(text=balanced)
        except Exception as e:
            messagebox.showerror("oops!", f"invalid equation :(\n{e}")
            self.output.config(text="")

class StoichiometryGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🧪 stoichiometry calculator")
        self.root.geometry("700x550")
        self.root.resizable(False, False)
        tk.Label(root, text="enter your unbalanced chemical equation:", font=("Consolas", 12)).pack(pady=10)
        eq_frame = tk.Frame(root)
        eq_frame.pack(pady=5)
    
        react_frame = tk.Frame(eq_frame)
        react_frame.grid(row=0, column=0, padx=10)
        tk.Label(react_frame, text="reactants", font=("Consolas", 10)).pack()
        self.react_entry = SubscriptEntry(react_frame, font=("Consolas", 12), width=25)
    
        self.react_entry.pack()
    
        tk.Label(eq_frame, text="→", font=("Consolas", 16)).grid(row=0, column=1, padx=5)
    
        prod_frame = tk.Frame(eq_frame)
        prod_frame.grid(row=0, column=2, padx=10)
        tk.Label(prod_frame, text="products", font=("Consolas", 10)).pack()
        self.prod_entry = SubscriptEntry(prod_frame, font=("Consolas", 12), width=25)
        self.prod_entry.pack()
    
        known_frame = tk.Frame(root)
        known_frame.pack(pady=15)
    
        tk.Label(known_frame, text="known compound:", font=("Consolas", 10)).grid(row=0, column=0, padx=5)
        self.known_compound_var = tk.StringVar()
        self.known_compound_menu = tk.OptionMenu(known_frame, self.known_compound_var, "")
        self.known_compound_menu.config(font=("Consolas", 10), width=15)
        self.known_compound_menu.grid(row=0, column=1, padx=5)
    
        amount_frame = tk.Frame(root)
        amount_frame.pack(pady=10)
    
        tk.Label(amount_frame, text="amount:", font=("Consolas", 10)).grid(row=0, column=0, padx=5)
        self.amount_entry = tk.Entry(amount_frame, font=("Consolas", 12), width=15)
        self.amount_entry.grid(row=0, column=1, padx=5)
    
        tk.Label(amount_frame, text="unit:", font=("Consolas", 10)).grid(row=0, column=2, padx=5)
        self.unit_var = tk.StringVar(value="grams")
        tk.OptionMenu(amount_frame, self.unit_var, "grams", "moles").grid(row=0, column=3, padx=5)
        unknown_frame = tk.Frame(root)
        unknown_frame.pack(pady=15)
    
        tk.Label(unknown_frame, text="find amount of:", font=("Consolas", 10)).grid(row=0, column=0, padx=5)
        self.unknown_compound_var = tk.StringVar()
        self.unknown_compound_menu = tk.OptionMenu(unknown_frame, self.unknown_compound_var, "")
        self.unknown_compound_menu.config(font=("Consolas", 10), width=15)
        self.unknown_compound_menu.grid(row=0, column=1, padx=5)
    
        tk.Label(unknown_frame, text="in:", font=("Consolas", 10)).grid(row=0, column=2, padx=5)
        self.output_unit_var = tk.StringVar(value="grams")
        tk.OptionMenu(unknown_frame, self.output_unit_var, "grams", "moles").grid(row=0, column=3, padx=5)
        btn_frame = tk.Frame(root)
        btn_frame.pack(pady=10)
    
        tk.Button(btn_frame, text="load equation", font=("Consolas", 11), bg="#FF9800", fg="white", command=self.load_equation).grid(row=0, column=0, padx=5)
    
        tk.Button(btn_frame, text="calculate", font=("Consolas", 11), bg="#2196F3", fg="white", command=self.calculate).grid(row=0, column=1, padx=5)
    
        self.output = tk.Label(root, text="", font=("Consolas", 12), fg="blue", wraplength=650, justify="left")
        self.output.pack(pady=20)
    
        self.balanced_eq = None
        self.coeffs = None
        self.compounds = None
    
    def load_equation(self):
        reactants = self.react_entry.get().strip()
        products = self.prod_entry.get().strip()
    
        if not reactants or not products:
            messagebox.showwarning("hey!", "you gotta enter both reactants and products!")
            return
    
        try:
            self.balanced_eq, self.coeffs, self.compounds = balance_equation(reactants, products)
            menu = self.known_compound_menu["menu"]
            menu.delete(0, "end")
    
            menu2 = self.unknown_compound_menu["menu"]
            menu2.delete(0, "end")
            for compound in self.compounds:
                menu.add_command(label=format_compound(compound), command=lambda c=compound: self.known_compound_var.set(c))
                menu2.add_command(label=format_compound(compound), command=lambda c=compound: self.unknown_compound_var.set(c))
    
            if self.compounds:
                self.known_compound_var.set(self.compounds[0])
                self.unknown_compound_var.set(self.compounds[0])
    
            self.output.config(text=f"balanced equation: {self.balanced_eq}")
    
        except Exception as e:
            messagebox.showerror("oops!", f"invalid equation :(\n{e}")
            self.output.config(text="")
    
    def calculate(self):
        if not self.balanced_eq:
            messagebox.showwarning("hey!", "load an equation first!")
            return
    
        known_compound = self.known_compound_var.get()
        unknown_compound = self.unknown_compound_var.get()
    
        if not known_compound or not unknown_compound:
            messagebox.showwarning("hey!", "select both compounds!")
            return
    
        try:
            amount = float(self.amount_entry.get())
        except ValueError:
            messagebox.showerror("oops!", "enter a valid number for amount!")
            return
    
        try:
            known_idx = self.compounds.index(known_compound)
            unknown_idx = self.compounds.index(unknown_compound)
    
            known_coeff = self.coeffs[known_idx]
            unknown_coeff = self.coeffs[unknown_idx]
    
            known_mm = calculate_molar_mass(known_compound)
            unknown_mm = calculate_molar_mass(unknown_compound)
            known_moles = amount / known_mm if self.unit_var.get() == "grams" else amount
            unknown_moles = known_moles * (unknown_coeff / known_coeff)
            result = unknown_moles * unknown_mm if self.output_unit_var.get() == "grams" else unknown_moles
            unit = "g" if self.output_unit_var.get() == "grams" else "mol"
            output_text = f"balanced equation: {self.balanced_eq}\n\n"
            output_text += f"given: {amount} {self.unit_var.get()} of {format_compound(known_compound)}\n"
            output_text += f"molar mass of {format_compound(known_compound)}: {known_mm:.2f} g/mol\n"
            output_text += f"molar mass of {format_compound(unknown_compound)}: {unknown_mm:.2f} g/mol\n\n"
            output_text += f"result: {result:.4f} {unit} of {format_compound(unknown_compound)}"
    
            self.output.config(text=output_text)
    
        except Exception as e:
            messagebox.showerror("oops!", f"calculation error :(\n{e}")

def main():
    root = tk.Tk()
    MainMenuGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()