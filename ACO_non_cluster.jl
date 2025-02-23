include("data_structures.jl")


module ACOInitialization
    
    import ..DataStructures: Instance, Patient, Nurse 
    using Random
    using Distributions

    struct AntColony
        pheromones::Matrix{Float64}
    end

    function initializeAntColony(instance::Instance, nurse_id::Int)::AntColony
        num_patients = length(instance.patients)
        pheromones = fill(1.0, num_patients, num_patients)  # Default pheromone values
        
        # Get nurse's starting location (Depot)
        depot_x, depot_y = instance.depot_x, instance.depot_y

        # Compute distance of each patient from depot
        patient_distances = [hypot(instance.patients[i].x_coord - depot_x, instance.patients[i].y_coord - depot_y) for i in 1:num_patients]
        
        # Normalize distances
        max_dist = maximum(patient_distances) + 1e-6
        normalized_distances = [1 - (d / max_dist) for d in patient_distances]  # Closer = higher pheromone

        # Assign initial pheromones based on proximity to depot
        for i in 1:num_patients
            for j in 1:num_patients
                pheromones[i, j] *= (normalized_distances[i] + normalized_distances[j]) / 2
            end
        end

        return AntColony(pheromones)
    end


    

    function updatePheromones!(colony::AntColony, instance::Instance, nurse_id::Int, selected_patient::Int, evaporation::Float64)
        # Increase pheromones for selected patient
        colony.pheromones[nurse_id, selected_patient] *= (1.0 - evaporation) + evaporation * 2.0

        # Reinforce pheromones for nearby patients
        for neighbor in 1:length(instance.patients)
            if neighbor != selected_patient
                dist = instance.travel_times[selected_patient, neighbor]
                if dist < 20  # Only reinforce close neighbors
                    colony.pheromones[nurse_id, neighbor] *= 1.05
                end
            end
        end
    end

    function computeSelectionProbabilities(colony::AntColony, instance::Instance, nurse_id::Int, possible_patients::Vector{Int}, alpha::Float64, beta::Float64)
        max_distance = maximum(instance.travel_times)
        probabilities = [
            (colony.pheromones[nurse_id, p] ^ alpha) *
            ((1.0 / (instance.travel_times[nurse_id, p] + 1e-6)) ^ beta) *
            exp(-instance.travel_times[nurse_id, p] / max_distance)  # Spatial penalty
            for p in possible_patients
        ]
        return probabilities / sum(probabilities)  # Normalize
    end

    function adaptiveEvaporation!(colony::AntColony, instance::Instance, nurse_id::Int)
        for p in 1:length(instance.patients)
            num_close_neighbors = sum(instance.travel_times[p, :] .< 20)  # Count nearby patients
            if num_close_neighbors > 5
                colony.pheromones[nurse_id, p] *= 0.95  # Faster evaporation in dense areas
            else
                colony.pheromones[nurse_id, p] *= 0.98  # Slower evaporation in sparse areas
            end
        end
    end

    function assignPatientsWithACO!(instance::Instance, alpha::Float64=1.0, beta::Float64=2.0, evaporation::Float64=0.5)
        num_nurses = length(instance.nurses)
        num_patients = length(instance.patients)

        # Step 1: Initialize pheromones with depot-based clustering
        colonies = [initializeAntColony(instance, nurse_id) for nurse_id in 1:num_nurses]
        
        nurse_loads = fill(0.0, num_nurses)
        unassigned_patients = Set(p.id for p in instance.patients)

        while !isempty(unassigned_patients)
            for (nurse_id, colony) in enumerate(colonies)
                if isempty(unassigned_patients)
                    break
                end
                
                # Get possible patients (respecting capacity)
                possible_patients = [p for p in unassigned_patients if nurse_loads[nurse_id] + instance.patients[p].demand <= instance.nurses[nurse_id].capacity]
                
                if isempty(possible_patients)
                    continue
                end
                
                # Compute selection probabilities
                probabilities = computeSelectionProbabilities(colony, instance, nurse_id, possible_patients, alpha, beta)
                
                # Select patient probabilistically
                selected_patient = possible_patients[rand(Categorical(probabilities))]
                push!(instance.nurses[nurse_id].route, selected_patient)
                nurse_loads[nurse_id] += instance.patients[selected_patient].demand
                delete!(unassigned_patients, selected_patient)

                # Update pheromones
                updatePheromones!(colony, instance, nurse_id, selected_patient, evaporation)
            end

            # Apply adaptive evaporation
            for (nurse_id, colony) in enumerate(colonies)
                adaptiveEvaporation!(colony, instance, nurse_id)
            end
        end
    end


end