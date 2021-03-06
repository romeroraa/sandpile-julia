# Sandpile64.jl
module Sandpile64

export sandpile_init64, run_sandpile64!, avalanche64!, add_grain!, is_unstable64

function sandpile_init64(x::Int64, y::Int64, setup::String="random")
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

function add_grain64!(z::Array{Int64,2})
    """
            add_grain!(z)
    Adds a grain to a random location on the  2d array grid for the Abelian sandpile model, `z`.
    Arguments:
            z (Array{Int64, 2}): 2D sandpile grid.
    """
    return z[rand(1:size(z)[1]), rand(1:size(z)[2])] += 1
end

function is_unstable64(z::Array{Int64,2}, fᵪ::Int64=4)
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

function avalanche64!(z::Array{Int64,2}, fᵪ=4)
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

function run_sandpile64!(z::Array{Int64,2}, N::Int64; N_crit::Bool=false, fᵪ::Int64=4)
    """
            run_sandpile(z, N, fᵪ)
    Run a sandpile model on an input grid `z` by adding `N` grains.
    Arguments:
            z (Array{Int64, 2}): 2D sandpile grid.
            N (Int64): Number of grains to add
            N_crit (boolean): Do you want to to all added `N` grains to result in collapses?
                              If set to `false`, `N` grains will be added to system but
                              returned sizes will only be those that result in avalanches.
                              Warning:  longer runtimes if set to `true`.
            fₓ (Int64): Critical value for sandpile model (Must be > 1)
    Returns
            (Vector{Int64}): Sizes of avalanches for each grain added.
    """
    s_list = zeros(Int64, N)
    if N_crit
            i = 1
            while i <= N
                add_grain64!(z)
                s = 0
                while is_unstable64(z, fᵪ)
                    s += avalanche64!(z)
                end
                if s > 0
                    s_list[i] = s
                    i += 1
                end
            end
    else
        for i in 1:N
            add_grain64!(z)
            s = 0
            while is_unstable64(z, fᵪ)
                s += avalanche64!(z)
            end
            s_list[i] = s
        end
    end
    return s_list[s_list .> 0]
end
end