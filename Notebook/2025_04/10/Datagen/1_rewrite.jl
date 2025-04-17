# save this as process_extxyz.jl and run with: julia process_extxyz.jl
cd("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_04/10/Datagen")
function process_extxyz(input_file::String, output_file::String)
    lines = readlines(input_file)
    i = 1
    new_frames = []

    while i <= length(lines)
        natoms = parse(Int, lines[i])
        comment = lines[i + 1]
        atoms = lines[i + 2 : i + 1 + natoms]

        push!(new_frames, (natoms, comment, atoms))
        i += 2 + natoms
    end

    # combine every two monomers into a dimer frame
    open(output_file, "w") do io
        for j in 1:2:length(new_frames)-1
            frame1 = new_frames[j]
            frame2 = new_frames[j+1]

            atoms1 = frame1[3]
            atoms2 = frame2[3]

            total_atoms = length(atoms1) + length(atoms2)
            println(io, total_atoms)
            println(io, "")  # no Lattice

            # Function to keep only first 4 columns: element, x, y, z
            clean_line = line -> join(split(line)[1:4], " ")

            for atom_line in atoms1
                println(io, clean_line(atom_line))
            end
            for atom_line in atoms2
                println(io, clean_line(atom_line))
            end
        end
    end
end
# Run it
process_extxyz("dimer_thiophene.extxyz", "dimer_thiophene.xyz")


dimer = read_extxyz("dimer_thiophene.extxyz")[1]
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

# Save as test1.xyz
process_xyz("test1.xyz", dimer)