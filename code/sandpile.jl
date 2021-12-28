# Sandpile.jl
# Author: Romero, Roland Albert
# These are a set of functions used to run basic Abelian sandpile models
# from Bak's Self-organized criticality (1988) paper
#
# Example Usage:
# x = 50
# y = 50
# z = sandpile_init(x, y, "zero");
# s = run_sandpile!(z, 100_000, 4);
#
#

function sandpile_init(x::Int64, y::Int64, setup::String="random")
    """
            sandpile_init(X, Y, setup)
    Initialize a 2d sandpile grid with size (`X`, `Y`). Can either have a zero or random
    initial condition.
    Arguments:
            x (Int64): length of grid along x
            y (Int64): length of grid along y
            setup (String): Type of inital setup for grid. Can either be `zero` or `random`.
    Returns:
            z (Matrix{Int64}): sandpile grid z with dimensions (x,y).
    """
    if setup == "random"
        return rand(0:3, x, y)
    elseif setup == "zero"
        return zeros(Int64, x, y)
    else
        error("`setup` must have value of `random` or `zero`")
    end
end

function add_grain!(z::Array{Int64,2})
    """
            add_grain!(z)
    Adds a grain to a random location on the  2d array grid for the Abelian sandpile model, `z`.
    Arguments:
            z (Array{Int64, 2}): 2D sandpile grid.
    """
    return z[rand(1:size(z)[1]), rand(1:size(z)[2])] += 1
end

function is_unstable(z::Array{Int64,2}, fᵪ::Int64=4)
    """
            is_unstable(z, fᵪ)
    Checks whether or not any location on the grid is unstable (has value above critical value).
    Arguments:
            z (Array{Int64, 2}): 2D sandpile grid.
            fₓ (Int64): Critical value for sandpile model (Must be > 1)
    Returns:
            (bool): Whether grid is stable or unstable.
    """
    if fᵪ <= 1
        error("`fᵪ` must be int greater than  1")
    end
    if length(findall(>=(fᵪ), z)) > 0
        return true
    else
        return false
    end
end

function avalanche!(z::Array{Int64,2}, fᵪ=4)
    """
            avalanche(z, fᵪ)
    Update grid to follow rules on locations that are unstable. 
    Follow Eq 3.2 of Self-organized criticality by Bak, Tang (1988).
    Returns how many avalanches (or slides s) occur for input state after collapsing all 
    unstable locations.
    Arguments:
            z (Array{Int64, 2}): 2D sandpile grid.
            fₓ (Int64): Critical value for sandpile model (Must be > 1)
    Returns:
            (Int64): Number of unstable locations collapsed.
    """
    s::Int64 = 0
    unstable_locs = findall(>=(fᵪ), z)
    for crit_loc in unstable_locs
        z[crit_loc] -= 4
        s += 1
        neighbors = [
            CartesianIndex(0, 1),
            CartesianIndex(0, -1),
            CartesianIndex(1, 0),
            CartesianIndex(-1, 0),
        ]
        for neighbor in neighbors
            try
                z[crit_loc + neighbor] += 1
            catch
            end
        end
    end
    return s
end

function run_sandpile!(z::Array{Int64,2}, N::Int64, fᵪ=4)
    """
            run_sandpile(z, N, fᵪ)
    Run a sandpile model on an input grid `z` by adding `N` grains.
    Arguments:
            z (Array{Int64, 2}): 2D sandpile grid.
            N (Int64): Number of grains to add
            fₓ (Int64): Critical value for sandpile model (Must be > 1)
    Returns
            (Vector{Int64}): Lifetimes of avalanches for each grain added.
    """
    s_list = zeros(Int64, N)
    for i in 1:N
        add_grain!(z)
        s = 0
        while is_unstable(z, fᵪ)
            s += avalanche!(z)
        end
        s_list[i] = s
    end
    return s_list
end
