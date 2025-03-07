
include("data_structures.jl")
include("ant_cluster.jl")
include("visualization.jl")
include("mutation.jl")
include("crossover.jl")

using Main.MyDataStructures

using .Visualization


function computeTotalTravelLength(instance::MyDataStructures.Instance)::Float64
    total_travel_length = 0.0
    for nurse in instance.nurses
        if length(nurse.route) > 0
            # Start from the depot
            prev_location = 0  # Assuming depot is index 0
            for patient_id in nurse.route
                total_travel_length += instance.travel_times[prev_location + 1, patient_id + 1]
                prev_location = patient_id
            end
            # Return to the depot
            total_travel_length += instance.travel_times[prev_location + 1, 1]
        end
    end
    return total_travel_length
end

"""
    computeFitness(instance::MyDataStructures.Instance) -> Float64

Evaluates the quality of an instance by summing over all nurses:
  - The total travel time along the route (as calculated by the ACO ordering).
  - Penalties for constraint violations:
      • A patient is served after its end_time.
      • The nurse returns to the depot after depot_return_time.
      • The nurse’s capacity is exceeded.
      
For each nurse, the function simulates the route:
  - Starting at time 0 at the depot (node 0).
  - For each patient in the nurse’s route, the travel time is added.
    If the arrival time is before the patient’s start_time, the nurse waits.
    If the arrival time is after the patient’s end_time, a penalty is incurred.
  - The nurse’s capacity usage is summed, and if it exceeds the nurse’s capacity, a penalty is added.
  - After the route, the travel time back to the depot is added, and if the return time exceeds depot_return_time, a penalty is incurred.
  
Finally, the fitness is defined as the inverse of (1 + total_cost), so that lower cost (fewer constraint violations and lower travel time) gives higher fitness.
"""
function calculateFitness(instance::MyDataStructures.Instance)
    # Penalty factors (adjust as needed)
    
    penalty_time = 100000.0   # penalty per time unit of lateness (either patient end_time or depot return)
    penalty_capacity = 1000.0  # penalty per unit of capacity exceeded
    
    total_travel_time = 0.0
    total_penalty = 0.0

    # Process each nurse's route
    for nurse in instance.nurses
        current_time = 0.0    # time starts at 0 at the depot
        current_node = 0      # depot is node 0 (corresponds to row/column 1 in travel_times)
        route_travel_time = 0.0
        capacity_used = 0.0

        for patient_id in nurse.route
            # Travel from current_node to the next patient
            travel = instance.travel_times[current_node + 1, patient_id + 1]
            route_travel_time += travel
            arrival_time = current_time + travel

            # Retrieve patient information
            patient = instance.patients[patient_id]

            # Determine if waiting is needed (waiting is allowed)
            waiting = 0.0
            if arrival_time < patient.start_time
                waiting = patient.start_time - arrival_time
                arrival_time = patient.start_time  # nurse starts care at patient.start_time
            end

            # Check if the arrival is too late (beyond the end_time)
            if arrival_time > patient.end_time
                total_penalty += (arrival_time - patient.end_time) * penalty_time
            end

            # Update current_time: add waiting and the care time
            current_time = arrival_time + patient.care_time
            capacity_used += patient.demand
            current_node = patient_id
        end

        # Add the travel time for the return trip to the depot
        return_travel = instance.travel_times[current_node + 1, 1]
        route_travel_time += return_travel
        arrival_return = current_time + return_travel

        # Check depot return constraint
        if arrival_return > instance.depot_return_time
            total_penalty += (arrival_return - instance.depot_return_time) * penalty_time
        end

        # Check nurse capacity constraint
        if capacity_used > nurse.capacity
            total_penalty += (capacity_used - nurse.capacity) * penalty_capacity
        end
        nurse.current_demand = capacity_used

        total_travel_time += route_travel_time
    end

    # Total cost is the travel time plus penalties.
    total_cost = total_travel_time + total_penalty

    # Define fitness so that lower total cost gives higher fitness.
    fitness = 1.0 / (1.0 + total_cost)
    instance.fitness = fitness
    return fitness
end

function applyRandomMutation!(instance::MyDataStructures.Instance)
    r = rand()
    if r < 0.33
        mutation_swap_one!(instance)
    elseif r < 0.85
        # Here, block_length is set to 2; you can tweak or randomize it.
        mutation_swap_multiple!(instance; block_length = 2)
    else
        mutation_move_one!(instance)
    end
end

# Select two random parents from the population.
function selectRandomParents(population::Vector{MyDataStructures.Instance})
    n = length(population)
    if n < 2
        error("Population must contain at least two individuals.")
    end
    idxs = randperm(n)[1:2]
    return population[idxs[1]], population[idxs[2]]
end

function runGA(initial_instance::MyDataStructures.Instance;
                           num_generations::Int=50, population_size::Int=10)
    # Generate the initial population by deep-copying the initial solution.
    depot_return_time = initial_instance.depot_return_time
    population = [deepcopy(initial_instance) for _ in 1:population_size]
    best_fitness = -Inf
    best_individual = nothing

    

    for generation in 1:num_generations
        fitnesses = [calculateFitness(instance) for instance in population]

        for (ind, fit) in zip(population, fitnesses)
            if fit > best_fitness
                best_fitness = fit
                best_individual = deepcopy(ind)
            end
        end

        sorted_indices = sortperm(fitnesses, rev=true)
    
        num_survivors = max(1, Int(round(population_size / 2)))
        survivors = [deepcopy(population[i]) for i in sorted_indices[1:num_survivors]]

        children = []
        for s in survivors
            parent1, parent2 = selectRandomParents(survivors)
            #child = crossoverInstances(parent1, parent2)
            child = deepcopy(parent1)
            applyRandomMutation!(child)
            
            push!(children, child)
        end

        # for child in children
        #     for 
        #     child, _ = constructRouteStrict(child.)
        # end

        # New population is survivors plus children.
        new_population = vcat(survivors, children)

        # Ensure population size is maintained.
        while length(new_population) < population_size
            push!(new_population, deepcopy(best_individual))
        end
        population = new_population[1:population_size]
        println("Generation $generation, best fitness = $best_fitness")
    end

    return best_individual, best_fitness
end


# Check the constraint violations for a single nurse.
function checkNurseConstraintViolations(nurse::MyDataStructures.Nurse, instance::MyDataStructures.Instance)
    violations = String[]
    current_time = 0.0
    current_node = 0  # Start at depot (node 0)
    capacity_used = 0.0

    for patient_id in nurse.route
        # Travel from current_node to next patient.
        travel_time = instance.travel_times[current_node + 1, patient_id + 1]
        current_time += travel_time

        # Retrieve patient info.
        patient = instance.patients[patient_id]

        # If arriving too early, wait until the patient's start time.
        if current_time < patient.start_time
            current_time = patient.start_time
        end

        # Check if nurse is too late to start patient care.
        if current_time > patient.end_time
            push!(violations, "Late arrival at patient $(patient.id) (arrival: $(current_time), deadline: $(patient.end_time))")
        end

        # Add the care time.
        current_time += patient.care_time

        # Update capacity usage and check capacity constraint.
        capacity_used += patient.demand
        if capacity_used > nurse.capacity
            push!(violations, "Capacity exceeded (used: $(capacity_used), capacity: $(nurse.capacity))")
        end

        # Set current node to this patient.
        current_node = patient_id
    end

    # Add travel time back to depot.
    travel_time = instance.travel_times[current_node + 1, 1]
    current_time += travel_time

    # Check depot return time constraint.
    if current_time > instance.depot_return_time
        push!(violations, "Depot return exceeded (arrival: $(current_time), depot return time: $(instance.depot_return_time))")
    end

    return violations
end

# Loop over all nurses in the instance and print constraint violations.
function checkAllNursesConstraints(instance::MyDataStructures.Instance)
    for nurse in instance.nurses
        violations = checkNurseConstraintViolations(nurse, instance)
        if !isempty(violations)
            println("Nurse $(nurse.id) violates constraints:")
            for v in violations
                println("  - ", v)
            end
        else
            println("Nurse $(nurse.id) meets all constraints.")
        end
    end
end



function main()
    instance = MyDataStructures.parseInstance("train/train_9.json")
    println("Parsed instance successfully!")

    assignPatientsWithACO!(instance)
    println("Assigned patients using Ant Colony Optimization:")
    for nurse in instance.nurses
        println("Nurse $(nurse.id): $(nurse.route)")
    end
    
    best_instance, best_fitness = runGA(instance, num_generations=100000, population_size=50)
    println("Best fitness: ", best_fitness)
    println("Best solution (nurse routes):")
    for nurse in best_instance.nurses
        println("Nurse $(nurse.id): $(nurse.route)")
    end
    Visualization.plotNurseRoutes(best_instance)
    println(computeTotalTravelLength(best_instance))
    print(checkAllNursesConstraints(instance))
end

main()