```
Transform extxyz to .xyz
```

function process_xyz(filename, atoms::JuLIP.Atoms)
    open(filename, "w") do io
        println(io, length(atoms))  # First line: number of atoms
        println(io, "")  # Second line: comment
        for (z, x) in zip(atoms.Z, atoms.X)
            element = chemical_symbol(z)  # convert atomic number to symbol, e.g., 6 -> C
            println(io, "$(element) $(x[1]) $(x[2]) $(x[3])")
        end
    end
end