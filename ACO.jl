module ACOInitialization

using Random
using Distributions
using Main.MyDataStructures

export assignPatientsWithACO!

"""
    computeTravelTimeForSolution(instance, routes)

Helper function that calculates the total travel time for a set of routes.
Each route is a vector of patient IDs. Assumes the depot is node 0.
"""
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

end  # module ACO



# module ACOStrictConstructive

#     using Random
#     using ..MyDataStructures

#     # export constructRouteStrict

#     """
#         constructRouteStrict(cluster::Vector{Int}, instance::Instance; depot_return_time, alpha, beta, rho, iterations)

#     Constructs a route (ordering of patient IDs) for a given cluster using a strict ACO:
#     - An ant starts at the depot (time 0) and iteratively selects a patient.
#     - A candidate patient is only considered if current_time + travel_time is within its time window.
#     - No waiting is allowed.
#     - Pheromone and a heuristic (inverse travel time) determine selection probabilities.

#     Returns a tuple (best_route, best_cost) where best_cost is the total travel time.
#     """
#     function constructRouteStrict(cluster::Vector{Int}, instance::MyDataStructures.Instance;
#                                 depot_return_time::Float64, alpha::Float64=1.0, beta::Float64=2.0,
#                                 rho::Float64=0.6, iterations::Int=10)
#         n = length(cluster)
#         # Create a local pheromone matrix.
#         # We map depot to index 1 and each patient in the cluster to indices 2..n+1.
#         pheromone = fill(1.0, n+1, n+1)
        
#         best_route = nothing
#         best_cost = Inf

#         for iter in 1:iterations
#             route = Int[]
#             current_node = 0   # 0 means depot
#             current_time = 0.0
#             remaining = Set(cluster)
#             feasible = true

#             while !isempty(remaining)
#                 candidates = Int[]
#                 scores = Float64[]
#                 for pid in remaining
#                     # Calculate travel time from current_node (depot=0, else patient id) to candidate.
#                     travel = instance.travel_times[current_node + 1, pid + 1]
#                     arrival_time = current_time + travel
#                     patient = instance.patients[pid]
#                     # Only consider candidate if arrival_time is within the time window
#                     if arrival_time >= patient.start_time && arrival_time <= patient.end_time
#                         push!(candidates, pid)
#                         # Use heuristic = 1/travel (the shorter the travel, the better)
#                         η = 1.0 / (travel + 1e-6)
#                         # Pheromone: map current_node to depot index if needed.
#                         from_idx = current_node == 0 ? 1 : findfirst(x -> x == current_node, cluster) + 1
#                         to_idx = findfirst(x -> x == pid, cluster) + 1
#                         τ = pheromone[from_idx, to_idx]
#                         push!(scores, (τ^alpha) * (η^beta))
#                     end
#                 end

#                 if isempty(candidates)
#                     feasible = false
#                     break
#                 end

#                 total_score = sum(scores)
#                 probabilities = total_score > 0 ? [s/total_score for s in scores] : fill(1/length(scores), length(scores))
#                 r = rand()
#                 cum = 0.0
#                 selected = candidates[1]
#                 for (j, p) in enumerate(probabilities)
#                     cum += p
#                     if r <= cum
#                         selected = candidates[j]
#                         break
#                     end
#                 end

#                 # Update route and state.
#                 push!(route, selected)
#                 remaining = setdiff(remaining, [selected])
#                 travel = instance.travel_times[current_node + 1, selected + 1]
#                 current_time += travel
#                 # Check feasibility: arriving time must be within patient window.
#                 patient = instance.patients[selected]
#                 if current_time < patient.start_time || current_time > patient.end_time
#                     feasible = false
#                     break
#                 end
#                 # Service the patient.
#                 current_time += patient.care_time
#                 current_node = selected

#                 # (Optional) Local pheromone update could be applied here.
#             end

#             # After constructing a route, add the return leg to depot.
#             if feasible
#                 total_cost = 0.0
#                 time_sim = 0.0
#                 prev = 0
#                 for pid in route
#                     travel = instance.travel_times[prev + 1, pid + 1]
#                     total_cost += travel
#                     time_sim += travel
#                     patient = instance.patients[pid]
#                     # Feasibility was ensured during construction.
#                     time_sim += patient.care_time
#                     prev = pid
#                 end
#                 total_cost += instance.travel_times[prev + 1, 1]
                
#                 # Also ensure the route allows return to depot by depot_return_time.
#                 if time_sim + instance.travel_times[prev + 1, 1] > depot_return_time
#                     feasible = false
#                 end

#                 if feasible && total_cost < best_cost
#                     best_cost = total_cost
#                     best_route = copy(route)
#                 end

#                 # Global pheromone update along this route.
#                 prev = 0
#                 for pid in route
#                     from_idx = prev == 0 ? 1 : findfirst(x -> x == prev, cluster) + 1
#                     to_idx = findfirst(x -> x == pid, cluster) + 1
#                     pheromone[from_idx, to_idx] = (1 - rho)*pheromone[from_idx, to_idx] + rho*(1.0/total_cost)
#                     prev = pid
#                 end
#             end
#         end

#         return best_route, best_cost
#     end

# end  # module ACOStrictConstructive


# module ACOLocalSearchStrict

#     using Random
#     using ..MyDataStructures
#     using ACOStrictConstructive

#     #export optimizeRouteStrict

#     """
#         optimizeRouteStrict(cluster::Vector{Int}, instance::Instance; depot_return_time, alpha, beta, rho, iterations, ls_iterations)

#     Uses a two-phase approach:
#     1. Constructs an initial route with the strict ACO (Variant 1).
#     2. Applies a simple 2‑opt local search to further refine the route, ensuring that any change keeps the arrival times within each patient’s time window.

#     Returns (best_route, best_cost).
#     """
#     function optimizeRouteStrict(cluster::Vector{Int}, instance::MyDataStructures.Instance;
#                                 depot_return_time::Float64, alpha::Float64=1.0, beta::Float64=2.0,
#                                 rho::Float64=0.6, iterations::Int=10, ls_iterations::Int=20)
#         # Phase 1: Construct initial route.
#         route, cost = ACOStrictConstructive.constructRouteStrict(cluster, instance;
#                                                                 depot_return_time=depot_return_time,
#                                                                 alpha=alpha, beta=beta, rho=rho, iterations=iterations)
#         if route == nothing
#             return nothing, Inf
#         end
#         best_route = copy(route)
#         best_cost = cost

#         # Define a helper that simulates a route and checks feasibility.
#         function simulateRoute(route)
#             time_sim = 0.0
#             total_travel = 0.0
#             prev = 0
#             for pid in route
#                 travel = instance.travel_times[prev + 1, pid + 1]
#                 arrival = time_sim + travel
#                 patient = instance.patients[pid]
#                 # Strict: arrival must be within the patient’s time window.
#                 if arrival < patient.start_time || arrival > patient.end_time
#                     return Inf, false
#                 end
#                 time_sim = arrival + patient.care_time
#                 total_travel += travel
#                 prev = pid
#             end
#             # Include return to depot.
#             total_travel += instance.travel_times[prev + 1, 1]
#             if time_sim + instance.travel_times[prev + 1, 1] > depot_return_time
#                 return Inf, false
#             end
#             return total_travel, true
#         end

#         # Phase 2: 2-opt local search.
#         for iter in 1:ls_iterations
#             improved = false
#             for i in 1:length(best_route)-1
#                 for j in i+1:length(best_route)
#                     new_route = copy(best_route)
#                     new_route[i:j] = reverse(new_route[i:j])
#                     new_cost, feasible = simulateRoute(new_route)
#                     if feasible && new_cost < best_cost
#                         best_cost = new_cost
#                         best_route = copy(new_route)
#                         improved = true
#                     end
#                 end
#             end
#             if !improved
#                 break
#             end
#         end

#         return best_route, best_cost
#     end

# end  # module ACOLocalSearchStrict
