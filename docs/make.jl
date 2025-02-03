using QFiND
using Documenter

DocMeta.setdocmeta!(QFiND, :DocTestSetup, :(using QFiND); recursive=true)

makedocs(;
    modules=[QFiND],
    authors="Hideaki Takahashi <hide.kyoto.1020@gmail.com> and contributors",
    sitename="QFiND.jl",
    format=Documenter.HTML(;
        canonical="https://htkhsh.github.io/QFiND.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/htkhsh/QFiND.jl",
    devbranch="main",
)
