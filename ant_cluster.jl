
using Random
using Distributions
using Main.MyDataStructures


function computeTravelTimeForSolution(instance::MyDataStructures.Instance, routes::Vector{Vector{Int}})
    total = 0.0
    for route in routes
        if !isempty(route)
            prev = 0  # start at depot
            for p in route
                total += instance.travel_times[prev + 1, p + 1]
                prev = p
            end
            # add cost to return to depot
            total += instance.travel_times[prev + 1, 1]
        end
    end
    return total
end

"""
    assignPatientsWithACO!(instance; iterations, alpha, beta, rho, Q)

ACO-based clustering:
- Runs for a number of iterations.
- In each iteration, each nurse builds its route from the depot using a probabilistic selection
  (roulette wheel) based on the score (pheromone^α * (1/travel_time)^β).
- After all routes are constructed, the total travel time is computed.
- Global pheromone update is applied: pheromones evaporate and then edges used in the current solution
  receive an additional deposit proportional to Q/total_travel_time.
- The best solution found over all iterations is then assigned to the instance.
"""
function assignPatientsWithACO!(instance::MyDataStructures.Instance; iterations::Int=100, alpha::Float64=1.0, beta::Float64=2.0, rho::Float64=0.6, Q::Float64=1.0)
    n_nodes = size(instance.travel_times, 1)  # (#patients + depot)
    n_nurses = length(instance.nurses)
    
    # Initialize pheromone matrix (for all edges)
    pheromone = fill(1.0, n_nodes, n_nodes)
    
    # Build a dictionary of all patients
    all_patients = Dict{Int, MyDataStructures.Patient}()
    for patient in instance.patients
        all_patients[patient.id] = patient
    end

    best_routes = [copy(nurse.route) for nurse in instance.nurses]
    best_total_time = Inf

    for iter in 1:iterations
        # For each iteration, reinitialize nurse routes and state
        local_routes = [Int[] for i in 1:n_nurses]
        current_nodes = fill(0, n_nurses)          # all nurses start at depot (node 0)
        remaining_capacities = [nurse.capacity for nurse in instance.nurses]
        unassigned = deepcopy(all_patients)         # reset unassigned patients
        active = trues(n_nurses)                    # flag for each nurse if further assignments are possible

        # Build routes until no feasible patient remains
        while !isempty(unassigned) && any(active)
            for i in 1:n_nurses
                if !active[i]
                    continue
                end

                current_index = current_nodes[i] + 1
                # Collect feasible candidates for nurse i
                feasible = []
                for (pid, patient) in unassigned
                    if patient.demand <= remaining_capacities[i]
                        push!(feasible, patient)
                    end
                end

                if isempty(feasible)
                    active[i] = false
                    continue
                end

                # Compute ACO scores for each candidate:
                scores = Float64[]
                for candidate in feasible
                    candidate_index = candidate.id + 1
                    travel_time = instance.travel_times[current_index, candidate_index]
                    # Avoid division by zero (if travel_time==0, use a very high heuristic value)
                    eta = travel_time > 0 ? 1.0 / travel_time : 1e6
                    tau = pheromone[current_index, candidate_index]
                    push!(scores, (tau^alpha) * (eta^beta))
                end

                # Normalize scores to form a probability distribution
                total_score = sum(scores)
                probabilities = total_score > 0 ? [s/total_score for s in scores] : fill(1/length(scores), length(scores))

                # Perform roulette wheel selection
                r = rand()
                cum_sum = 0.0
                selected_index = 1
                for j in 1:length(probabilities)
                    cum_sum += probabilities[j]
                    if r <= cum_sum
                        selected_index = j
                        break
                    end
                end

                best_candidate = feasible[selected_index]

                # Update nurse i’s route and state
                push!(local_routes[i], best_candidate.id)
                remaining_capacities[i] -= best_candidate.demand
                old_node = current_nodes[i]
                current_nodes[i] = best_candidate.id
                delete!(unassigned, best_candidate.id)

                # Local pheromone update on the edge (optional)
                candidate_index = best_candidate.id + 1
                pheromone[current_index, candidate_index] = (1 - rho) * pheromone[current_index, candidate_index] +
                                                            rho * (1.0 / instance.travel_times[current_index, candidate_index])
                pheromone[candidate_index, current_index] = (1 - rho) * pheromone[candidate_index, current_index] +
                                                            rho * (1.0 / instance.travel_times[current_index, candidate_index])
            end
        end  # end while

        # Evaluate the solution built in this iteration
        total_time = computeTravelTimeForSolution(instance, local_routes)
        if total_time < best_total_time
            best_total_time = total_time
            best_routes = deepcopy(local_routes)
        end

        # Global pheromone update: evaporate pheromones on all edges and reinforce edges used in this iteration’s routes.
        for i in 1:n_nurses
            prev_node = 0
            for pid in local_routes[i]
                from_idx = prev_node + 1
                to_idx = pid + 1
                pheromone[from_idx, to_idx] = (1 - rho) * pheromone[from_idx, to_idx] + rho * (Q / total_time)
                pheromone[to_idx, from_idx] = (1 - rho) * pheromone[to_idx, from_idx] + rho * (Q / total_time)
                prev_node = pid
            end
            # Account for the return to the depot.
            from_idx = prev_node + 1
            to_idx = 1
            pheromone[from_idx, to_idx] = (1 - rho) * pheromone[from_idx, to_idx] + rho * (Q / total_time)
            pheromone[to_idx, from_idx] = (1 - rho) * pheromone[to_idx, from_idx] + rho * (Q / total_time)
        end
    end  # end iterations loop

    # Assign the best routes found to the instance.
    for i in 1:n_nurses
        instance.nurses[i].route = best_routes[i]
    end
    return
end