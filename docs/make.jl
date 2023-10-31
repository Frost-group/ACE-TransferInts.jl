using ACE-TransferInts
using Documenter

DocMeta.setdocmeta!(ACE-TransferInts, :DocTestSetup, :(using ACE-TransferInts); recursive=true)

makedocs(;
    modules=[ACE-TransferInts],
    authors="Jarvist Moore Frost",
    repo="https://github.com/Frost-group/ACE-TransferInts.jl/blob/{commit}{path}#{line}",
    sitename="ACE-TransferInts.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Frost-group.github.io/ACE-TransferInts.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Frost-group/ACE-TransferInts.jl",
    devbranch="main",
)
