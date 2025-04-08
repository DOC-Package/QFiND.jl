using QFiND
using Documenter

DocMeta.setdocmeta!(QFiND, :DocTestSetup, :(using QFiND); recursive=true)

makedocs(;
    modules=[QFiND],
    authors="Hideaki Takahashi <takahashi.hideaki.w33@kyoto-u.jp> and contributors",
    sitename="QFiND.jl",
    format=Documenter.HTML(;
        canonical="https://doc-package.github.io/QFiND.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Basic Theory" => "theory.md",
        "ID Discretization" => "id.md",
        "BSDO Discretization" => "bsdo.md",
        "Work Flow" => "workflow.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/DOC-Package/QFiND.jl",
    devbranch="main",
)
