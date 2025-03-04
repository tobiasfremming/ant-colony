

using Random

# Repair function for unscheduled patients.
# Ensures every unscheduled patient is assigned to exactly one nurse route.
function repairUnscheduled!(instance::MyDataStructures.Instance, unscheduled::Vector{Int})
    for p in unscheduled
        # Double-check that p is not already assigned.
        already_assigned = any(p in nurse.route for nurse in instance.nurses)
        if already_assigned
            continue
        end

        assigned = false

        # First, try to assign the patient to a nurse with an empty route.
        for nurse in instance.nurses
            if isempty(nurse.route)
                push!(nurse.route, p)
                println("Assigned unscheduled patient $p to nurse $(nurse.id) (empty route).")
                assigned = true
                break
            end
        end

        if assigned
            continue
        end

        # If no empty route is available, try to insert into existing routes.
        best_cost_increase = Inf
        best_insertion = nothing  # Will hold (nurse, insertion_position)
        for nurse in instance.nurses
            current_route = nurse.route
            # Consider every possible insertion position (from 1 to end+1).
            for pos in 1:(length(current_route) + 1)
                new_route = copy(current_route)
                insert!(new_route, pos, p)
                cost_new = route_cost(new_route, instance)
                cost_current = route_cost(current_route, instance)
                cost_increase = cost_new - cost_current
                if cost_increase < best_cost_increase
                    best_cost_increase = cost_increase
                    best_insertion = (nurse, pos)
                end
            end
        end

        if best_insertion !== nothing
            (selected_nurse, pos) = best_insertion
            insert!(selected_nurse.route, pos, p)
            println("Inserted unscheduled patient $p into nurse $(selected_nurse.id)'s route at position $pos.")
            assigned = true
        end

        # If, for some reason, no insertion was found, force assign to the first nurse.
        if !assigned
            println("Forcing insertion of unscheduled patient $p into nurse $(instance.nurses[1].id)'s route at the end.")
            push!(instance.nurses[1].route, p)
        end
    end
    return
end



# Mutation operator 1: Swap one patient between two nurses.
function mutation_swap_one!(instance::MyDataStructures.Instance)
    n_nurses = length(instance.nurses)
    if n_nurses < 2
        return
    end
    # Choose two distinct nurses randomly.
    idxs = randperm(n_nurses)[1:2]
    nurseA = instance.nurses[idxs[1]]
    nurseB = instance.nurses[idxs[2]]
    # Ensure both nurses have at least one patient.
    if isempty(nurseA.route) || isempty(nurseB.route)
        return  # Nothing to swap if one route is empty.
    end
    # Select a random patient from each nurse's cluster.
    posA = rand(1:length(nurseA.route))
    posB = rand(1:length(nurseB.route))
    # Swap the patients.
    nurseA.route[posA], nurseB.route[posB] = nurseB.route[posB], nurseA.route[posA]
end

# Mutation operator 2: Swap multiple patients between nurses.
# This performs a specified number of independent swaps.
# Mutation operator 2: Swap a contiguous block of patients between two nurses.
function mutation_swap_multiple!(instance::MyDataStructures.Instance; block_length::Int=2)
    n_nurses = length(instance.nurses)
    if n_nurses < 2
        return
    end
    # Select two distinct nurses at random.
    idxs = randperm(n_nurses)[1:2]
    nurseA = instance.nurses[idxs[1]]
    nurseB = instance.nurses[idxs[2]]
    
    # Ensure both nurses have enough patients for a block swap.
    if length(nurseA.route) < block_length || length(nurseB.route) < block_length
        return  # not enough elements to swap a block
    end
    
    # Choose random starting indices such that a block of the desired length fits.
    startA = rand(1:(length(nurseA.route) - block_length + 1))
    startB = rand(1:(length(nurseB.route) - block_length + 1))
    
    # Swap the blocks.
    blockA = nurseA.route[startA : startA + block_length - 1]
    blockB = nurseB.route[startB : startB + block_length - 1]
    for i in 0:(block_length - 1)
        nurseA.route[startA + i] = blockB[i + 1]
        nurseB.route[startB + i] = blockA[i + 1]
    end

end


# Mutation operator 3: Move one patient from one nurse's cluster to another.
function mutation_move_one!(instance::MyDataStructures.Instance)
    n_nurses = length(instance.nurses)
    if n_nurses < 2
        return
    end
    # Choose a donor nurse that has at least one patient.
    donor_index = rand(1:n_nurses)
    while isempty(instance.nurses[donor_index].route)
        donor_index = rand(1:n_nurses)
    end
    # Choose a recipient nurse (different from donor).
    recipient_index = rand(1:n_nurses)
    while recipient_index == donor_index
        recipient_index = rand(1:n_nurses)
    end
    donor = instance.nurses[donor_index]
    recipient = instance.nurses[recipient_index]
    # Remove a random patient from the donor's cluster.
    pos = rand(1:length(donor.route))
    patient = donor.route[pos]
    deleteat!(donor.route, pos)
    # Insert the patient into the recipient's cluster at a random position.
    pos_recipient = isempty(recipient.route) ? 1 : rand(1:length(recipient.route)+1)
    insert!(recipient.route, pos_recipient, patient)
end




# Add heurisitc mutation