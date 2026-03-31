from ase.lattice.cubic import FaceCenteredCubic
import numpy as np

# Criar uma estrutura base (ex: FCC - Comum em HEAs)
# Vamos criar um cubo 10x10x10 (4000 átomos no FCC)
elements = ['Co', 'Cr', 'Fe', 'Ni', 'Mn', 'Cu'] # Exemplo de UHEA
size = (10, 10, 10)
atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          size=size, symbol='Cu') # Começa tudo como Cu

# Distribuir os elementos aleatoriamente
n_atoms = len(atoms)
new_symbols = np.random.choice(elements, size=n_atoms)
atoms.set_chemical_symbols(new_symbols)

# Agora 'atoms' é o bloco de cristal aleatório.
