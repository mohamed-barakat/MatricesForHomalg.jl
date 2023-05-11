using MatricesForHomalg
using Documenter

DocMeta.setdocmeta!(MatricesForHomalg, :DocTestSetup, :(using MatricesForHomalg); recursive=true)

makedocs(;
    modules=[MatricesForHomalg],
    authors="Mohamed Barakat <mohamed.barakat@uni-siegen.de>, Johanna Knecht <johanna.knecht@student.uni-siegen.de>",
    repo="https://github.com/homalg-project/MatricesForHomalg.jl/blob/{commit}{path}#{line}",
    sitename="MatricesForHomalg.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://homalg-project.github.io/MatricesForHomalg.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/homalg-project/MatricesForHomalg.jl",
    devbranch="main",
)
