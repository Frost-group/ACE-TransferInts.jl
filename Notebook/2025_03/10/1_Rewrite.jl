# save this as process_extxyz.jl and run with: julia process_extxyz.jl

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
i = 18
process_extxyz("C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_03/10/test"*string(i)*".extxyz",
 "C:/Users/Keerati/Documents/GitHub/ACE-TransferInts.jl/Notebook/2025_03/10/test"*string(i)*".xyz")