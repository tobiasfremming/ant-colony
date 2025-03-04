
# Order Crossover (OX): Given two parent orderings (each a permutation of all patient IDs),
# produce a child ordering that preserves the relative order of patients.
function orderCrossover(p1::Vector{Int}, p2::Vector{Int})
    M = length(p1)
    child = fill(-1, M)
    # Choose two random cut points
    c1 = rand(1:M)
    c2 = rand(1:M)
    if c1 > c2
        c1, c2 = c2, c1
    end
    # Copy the segment from p1 to the child.
    for i in c1:c2
        child[i] = p1[i]
    end
    # Fill the remaining positions with the order of elements in p2 that are not yet in child.
    pos = mod(c2, M) + 1  # next position (wrap-around)
    for p in p2
        if !(p in child)
            child[pos] = p
            pos = mod(pos, M) + 1
        end
    end
    return child
end



# Crossover operator for two instances.
# It creates a new instance by combining the routes (patient orderings) of the two parents.
# The child is built by:
#  1. Flattening each parent's nurse routes into a single ordering.
#  2. Using order crossover to produce a new ordering.
#  3. Partitioning the new ordering into nurse routes (here we split evenly).
function crossoverInstances(parent1::MyDataStructures.Instance, parent2::MyDataStructures.Instance)
    # Create a child instance by copying parent1 (assumes both parents share same instance parameters).
    child = deepcopy(parent1)
    
    # Flatten the nurse routes into a global ordering.
    ordering1 = Int[]
    ordering2 = Int[]
    for nurse in parent1.nurses
        append!(ordering1, nurse.route)
    end
    for nurse in parent2.nurses
        append!(ordering2, nurse.route)
    end
    
    # Apply order crossover to produce a new permutation of patient IDs.
    child_order = orderCrossover(ordering1, ordering2)
    
    # Partition the child_order into nurse routes.
    n_nurses = length(child.nurses)
    M = length(child_order)
    # Evenly split the global ordering.
    base = div(M, n_nurses)
    rem = M % n_nurses  # remainder to distribute among first few nurses.
    index = 1
    for i in 1:n_nurses
        # Give nurse i an extra patient if needed.
        size = base + (i <= rem ? 1 : 0)
        # If a nurse gets an empty route (could happen if there are fewer patients than nurses),
        # assign an empty route.
        if size > 0
            child.nurses[i].route = child_order[index : index + size - 1]
        else
            child.nurses[i].route = Int[]
        end
        index += size
    end

    return child
end