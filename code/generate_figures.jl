# This script uses the functions in sandpile.jl
include("Sandpile.jl")
using .Sandpile
using Plots

x, y = 50, 50; # Dimensions of grid
N = 100_000; # Number of grains to be added
fᵪ = 4; # Critical value for sandpile model

for samples in [1, 10, 200]
    for init in ["zero", "random"]

        s = [run_sandpile!(sandpile_init(x, y, "zero"), N) for i in 1:samples]

        # Get counts of each slide size s
        unique_s = sort(unique(reduce(vcat, s)))
        counts = [count(==(element), reduce(vcat, s)) for element in unique_s]
        yticks = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
        xticks = [1e0, 1e1, 1e2, 1e3]
        plot(
            unique_s,
            counts ./ sum(counts);
            xaxis=:log,
            yaxis=:log,
            leg=false,
            yticks=yticks,
            xticks=xticks,
            size=(600, 600),
            dpi=300,
        )
        plot!(
            [1e0, 1e5],
            [1.25e-1, 1.25e-6];
            xaxis=:log,
            yaxis=:log,
            xlabel="s",
            ylabel="D(s)",
            xlim=(1e0, 1.5e3),
            ylim=(1e-5, 1.5e-1),
            linestyle=:dash,
            dpi=300,
        )
        savefig("../results/sandpile_dim$(x)x$(y)_N$(N)_init$(init)_$(samples)samples.png")

    end
end
