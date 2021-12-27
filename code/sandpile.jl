radius = 2
z = zeros(Int64, 2radius+1, 2radius+1)

function add_grain!(z::Array{Int64, 2})
        z[rand(1:size(z)[1]), rand(1:size(z)[2])] += 1
end


function is_unstable(z::Array{Int64, 2}; fᵪ=4)
        unstable_locs = findall(>=(fᵪ), z)
        if length(unstable_locs) > 0
                return true
        else
                return false
        end
end


function avalanche!(z::Array{Int64,  2}; fᵪ=4)
        lifetime = 0
        unstable_locs = findall(>=(fᵪ), z)
        for crit_loc in unstable_locs
                z[crit_loc] -= 4
                lifetime+=1
                neighbors = [CartesianIndex(0,1), CartesianIndex(0,-1),
                            CartesianIndex(1,0), CartesianIndex(-1,0)]
                for neighbor in neighbors
                        try
                                z[crit_loc+neighbor]+=1
                        catch
                        end
                end
        end
        return lifetime
end



fᵪ = 4
z = zeros(Int64, 50, 50)
lifetimes = zeros(Int64, 1000000)
for i in 1:length(lifetimes)
        add_grain!(z)
        lifetime = 0
        while is_unstable(z)
                lifetime += avalanche!(z)
        end
        lifetimes[i]=lifetime
end

using Plots
heatmap(1:50, 1:50, z)
savefig("test.png")
